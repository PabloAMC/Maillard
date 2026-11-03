# Zhou, Wang, Li, Wang, Sun, Chen, Dong & Ruan 2020 — EXTRACTION (file stem `zhou2020`)
DOI: 10.1002/jsfa.10327 (printed on p. 2982, top left; Crossref 2026-11-02: "Application of Maillard reaction product
of xylose–pea protein enzymatic hydrolysate in 3D printing", J. Sci. Food Agric. **100**(7):2982-2990, issued 2020-05,
authors Quancheng Zhou, Meng Wang, Hui Li, Sai Wang, Wanting Sun, Xuanhong Chen, Leichao Dong, Zheng Ruan. **Matches the
PDF** (title, all eight authors, journal, volume, pages; the PDF's metadata subject line reads "2020.100:2982-2990").
**The file is the JSFA paper the brief names.** Received 29 Oct 2019, accepted 12 Feb 2020, online 2 Mar 2020.)
**Not to be confused with** the Jiangnan / Ho dossiers `zhou2022`, `zhou2023`, `zhou2023c`, `zhou2024`, `zhou2024b`,
`zhou2025b` (first author **Tong Zhou**, Ala-xylose Amadori pots), nor with `Zhou2025`, `zhou2000`, `zhou2002`,
`zhou2026` (other Zhous). This first author is **Quancheng Zhou** (Shandong Univ. of Technology / Nanchang); no author
is shared with any dossier on disk.

### Pea protein hydrolysed with "neutral protease" (40.3 °C, pH 7.02, 3 h, 6.80 % protein), the clarified hydrolysate (50 mL) heated with 0 / 2 / 4 / 6 / 8 / 10 g xylose for 48 h in a 60 °C water bath, then diluted 20× and kneaded into a potato-starch / pea-protein / glycerol dough for extrusion 3D printing. **Measured: SEM, FTIR transmittance, XRD, texture (cohesion, adhesiveness), DSC Tg, and O2•− / •OH scavenging of the dough materials.** **No Maillard chemistry is measured on the reaction product itself: no A420, no A294, no colour coordinates, no pH, no free amino groups, no sugar loss, no volatiles, no time or temperature series.**

**Read 2026-11-02 from `data/articles/zhou2020.pdf`** (9 pp., Wiley/SCI, born-digital). Prose for all pages from
`pdftotext -layout` (scratchpad `ingest/zhou2020_work/zhou2020.txt`). All nine pages rendered at 150 dpi
(`ingest/zhou2020_work/p-1.png` … `p-9.png`); **pp. 2982 (methods), 2987 (Table 3, Fig. 5), 2988 (Fig. 6) and 2989
(Table 4, Figs. 7-8) read by eye from the renders**; **Fig. 8 digitised from a 300 dpi crop** (`fig8-8.png`, read error
±0.5 percentage points). Tables 1-2 (FTIR) and Figs. 1-4 (SEM, FTIR, XRD) seen in the text layer only; they are
functional / structural and not transcribed. pdfimages confirms the figures are raster (300 ppi).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure with
stated read error · **[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]**
derived by me, never printed · **[!]** flag.

**Laboratory.** Corresponding author **Zheng Ruan** (State Key Laboratory of Food Science and Technology, Institute of
Nutrition and School of Food Science, **Nanchang University**, Jiangxi; ruanzheng@ncu.edu.cn); first author
**Quancheng Zhou** holds both affiliations (Nanchang and **School of Agricultural Engineering and Food Science, Shandong
University of Technology, Zibo**); the other six authors are Shandong Univ. of Technology only. Funding: Open Project
Program SKLF-KF201817 of the Nanchang State Key Lab. The bench work is most likely Zhou QC's SDUT group (their own
earlier potato-starch / pea-protein 3D-printing papers are refs 2 and 12) [Z, ~75 %].
**Not in `data/keys/laboratories.yml`: NEW.** Proposed id **`nanchang_ruan_z`** (corresponding author first, per the
file's rule), lab string **"Nanchang Univ. / Shandong Univ. of Technology (Zhou QC / Ruan Z)"**. Not the Jiangnan
Zhang X / Ho network (`jiangnan_zhang_x`; that Zhou is Tong Zhou), and not the Jiangnan pea / xylose work. The surname
Ruan in `ruan2014` (Qijun Ruan, Jiangnan, Hua Y) and `ruan2025` (Henan Univ. of Technology, Yin W-t) is different
people. **The id is only needed if a ledger entry is made; sec. 7 proposes none.**

## 1. What this paper answers for the model

| model need | what Zhou 2020 gives | tag |
|---|---|---|
| amine source composition (free amino N, degree of hydrolysis, peptide MW, amino-acid profile) | **none**: enzyme dose, T, pH, time, protein % only; supernatant composition not measured | [M] / [!] |
| reducing-sugar charge | 2 / 4 / 6 / 8 / 10 g xylose per 50 mL hydrolysate ≈ 0.27-1.33 mol/L nominal | [M] / [Z] |
| T / t of the Maillard step | 60 °C, 48 h, water bath; one condition, no series | [M] |
| pH initial / final | hydrolysis set at pH 7.02; **no pH at xylose addition, during or after the 48 h** | [M] / [!] |
| browning (A420), intermediates (A294), colour (L*a*b*) | **none**; only photographs of the printed dough (Fig. 5) after 20× dilution into starch dough | [!] |
| amine loss (OPA / TNBS), sugar loss | **none** | [!] |
| volatiles (any lane) | **none** | [!] |
| cysteine / sulfur | **none added, none measured** | [M] |
| a number that depends on the Maillard extent | O2•− and •OH scavenging (%) of the dough materials vs xylose dose (Fig. 8): not a model observable | [F] |
| binding / partition (matrix layer) | none | — |
| functional (SEM, FTIR T %, XRD, TPA, DSC Tg) | Tables 1-4, Figs. 1-7; skipped (one-line notes in sec. 4) | [M] |

**Answers.**
1. **Exact Maillard conditions** (confidence ~90 % on what is printed; the gaps are real gaps, not misreadings):
   - **Hydrolysate**: commercial pea protein (Yantai Shuangta Food Co.; grade and purity not stated), 6.80 % protein
     (w/v assumed), **"neutral protease, 10.6 %"** (basis not stated; most likely enzyme : substrate w/w [Z, ~70 %]),
     40.3 °C, pH 7.02, 3 h, water bath; enzyme inactivated 10 min in boiling water; centrifuged, supernatant vacuum
     filtered [M]. The odd decimals (40.3 °C, 7.02, 6.80 %, 10.6 %) look like response-surface optima from an
     unreported optimisation [Z]. **No degree of hydrolysis, no MW distribution, no soluble protein / peptide content,
     no free amino N, no amino-acid profile** [!].
   - **Xylose**: 2, 4, 6, 8, 10 g mixed with 50 mL hydrolysate; CK = hydrolysate without xylose [M]. Nominal
     40-200 g/L = **0.27 / 0.53 / 0.80 / 1.07 / 1.33 mol/L** if 50 mL is the final volume; 2-11 % lower (down to
     ~1.19 mol/L for 10 g) if 50 mL is the volume before the sugar dissolves (xylose ~95 mL/mol) [Z].
   - **Amine side (unmeasured)**: at 68 g/L protein, if all of it were soluble, Lys ε-NH2 ~30 mM plus new α-NH2 from
     DH 5-20 % ~30-120 mM, i.e. **~0.03-0.15 mol/L total amine, xylose : amine ≈ 2-40 molar** [Z, order of magnitude
     only; protein purity, solubility and DH unknown]. Sugar is in large excess in every arm.
   - **Additives**: none (no cysteine, thiamine, buffer, salt) [M].
   - **pH**: 7.02 is the hydrolysis set-point; not re-adjusted, not buffered, not measured afterwards [M / !]. An
     unbuffered pentose pot at 60 °C for 48 h would be expected to acidify, but this paper gives no number.
   - **T / t / vessel / water**: 60 °C water bath, 48 h, "filtered and stored at 4 °C" [M]; vessel, closure and
     headspace **not stated**; fully aqueous (aw ≈ 1 by construction [Z]).
2. **Measured chemistry: none on the reaction product** (confidence ~97 %). No A420, A294, colour meter, pH, free
   amino groups, xylose, Amadori, furfural or any volatile appears in methods, tables or figures (`grep` of the text
   layer for A420 / absorbance / brown / pH / colour / amino returns only the hydrolysis set-point and the antioxidant
   assays' absorbance readings). Colour is discussed only qualitatively from photographs (Fig. 5) [Q]; "chromaticity"
   is named as an index for *future* work [Q]. The only Maillard-dependent numbers are the dough's radical scavenging
   rates against the xylose dose (Fig. 8) and structural readouts of the dough (FTIR, XRD, texture, DSC). **Nothing
   kinetic; a single 60 °C / 48 h endpoint per dose, and even that endpoint is not characterised chemically.**
3. **Usable as a held-out "pea pot"? No.** P(usable for the trunk lane) ~1 %; P(usable for the sulfur lane) ~0 % (no
   sulfur source, no sulfur readout). This is a 3D-printing functionality paper whose only "MRP characterisation" is
   indirect, measured on a starch-dominated dough after 20× dilution, drying at 38 °C and grinding (sec. 2). Even
   qualitatively it gives the model only one thing: a stated recipe (pea hydrolysate + 0.27-1.33 M xylose, 60 °C,
   48 h) that a designed experiment could reuse [P]. It is not a pea Maillard data set.
4. **Proposed ledger entries: none** (confidence ~95 % that this is right; sec. 7 gives the reasons and the
   precedent). Fig. 8 is quantitative but measures antioxidant capacity of a printed dough, an observable no lane
   predicts.

## 2. Conditions

- **Proteins / reactants** [M]: pea protein (Yantai Shuangta Food Co. Ltd; isolate vs concentrate, protein %, and lot
  not stated); xylose, neutral protease and edible glycerol (Shanghai Ichun Biotechnology Co. Ltd; purity and enzyme
  activity not stated); potato starch "ordinary commercial material".
- **Hydrolysis** [M]: neutral protease 10.6 %; 40.3 °C; pH 7.02 (adjusting agent not stated); pea protein 6.80 %;
  3 h, water bath; inactivation 10 min in boiling water; centrifugation (speed / time not stated); supernatant vacuum
  filtered = "enzyme hydrolysate". [!] The boiling step is a 10 min ~100 °C heating of the hydrolysate **without**
  xylose (CK and all arms share it).
- **Maillard step** [M]: xylose 2 / 4 / 6 / 8 / 10 g + 50 mL hydrolysate, "reacted for 48 h in a 60 °C water bath"
  [V]; filtered; stored at 4 °C. Arms named "2 g xMPRs" … "10 g xMPRs" (sic); CK = no xylose. **No pH before / after,
  no buffer, no vessel, no closure, no headspace / O2 control, no heat-up, no replicate count.**
- **Dough** [M]: potato starch + pea protein at 1 % of starch mass, blended; xMRPs "diluted with distilled water 20
  times by weight were added" [V]; gelatinised at 67 °C (time not stated); 2 g edible glycerol added after
  gelatinisation; kneaded. [!] **The amount of diluted xMRP added, and the starch mass, are not stated**, so the xMRP
  fraction of the dough cannot be computed.
- **Printing / post-treatment** [M]: PORIMY printer (Wuxi), 1.0 mm circular nozzle, 15 mm/s, layer height 1.0 mm,
  heart-with-keyhole shape; printing temperature not stated. Materials "dried before and after printing at 38 °C" to
  constant weight, pulverised, 100-mesh sieved, stored in sealed bags.
- **Replicates / statistics** [M]: mean ± SD, one-way ANOVA + Duncan (SPSS 19), t-test in Excel, P < 0.05; **n not
  stated** [!].

## 3. Analytical methods

- **Antioxidant assays** [M] (the only solution chemistry): •OH scavenging by the o-phenanthroline / Fenton method
  (printed reagent names "o-diazophane" and "ferrous sulfide" are almost certainly o-phenanthroline and ferrous sulfate
  [Z, ~90 %]), 0.02 M phosphate pH 7.4, H2O2 0.02 %, 37 °C / 60 min, A536; O2•− scavenging by pyrogallol
  autoxidation, Tris-HCl pH 8.2, 37 °C, A320 at 30 s intervals. Scavenging % from the printed formulas. [!] The
  methods specify "xMRP solution" as the sample, but Fig. 8 is titled and discussed as the activity of the printing
  **materials** before and after printing; which was assayed (and at what dilution or extraction) is unresolved.
- **Structural** [M]: SEM (FEI Quanta 250, 10 kV, ×5000); FTIR (KBr 1:100, 400-4000 cm⁻¹); XRD (Bruker D8-Advance,
  3-50°, 0.02° steps); TPA (70 % compression, 1.0 mm/s test speed, 5 s interval); DSC (6-7 mg, sealed Al pan,
  10 °C/min, 0-180 °C, N2 50 mL/min).
- **No quantification of any Maillard species or index**; no standards, LODs or internal standards apply.
- **What the ± means**: SD, n unstated (Fig. 6, 8 error bars). Tables 1-4 print no ± at all.

## 4. Data

### 4a. Model-relevant data: none

No table or figure carries a Maillard chemistry readout (sec. 1, answer 2).

### 4b. Fig. 8 (p. 2989): radical scavenging of the dough materials vs xylose dose [F], ±0.5 points, 300 dpi crop

Recorded for completeness: it is the only xylose-dose series in the paper with a chemical (not mechanical) readout,
but it is an antioxidant-capacity index of a diluted, dried starch dough, not a model observable. x-axis printed
"xylose (/g)"; 0 = CK.

| xylose (g / 50 mL) | 0 | 2 | 4 | 6 | 8 | 10 |
|---|---|---|---|---|---|---|
| O2•− scavenging %, before printing | 28.5 | 33.6 | 37.5 | 42.0 | 40.4 | 38.3 |
| O2•− scavenging %, after printing | 33.3 | 37.2 | 40.9 | 45.3 | 44.3 | 42.9 |
| •OH scavenging %, before printing | 30.2 | 39.3 | 42.4 | 45.6 | 45.0 | 42.6 |
| •OH scavenging %, after printing | 35.6 | 43.4 | 47.0 | 49.6 | 49.1 | 47.1 |

Error bars ≈ ±1 point throughout [F]. Maximum at 6 g in all four series, as the caption states [M/Q]. The fall above
6 g cannot be interpreted chemically: the dough's xMRP loading is unknown and the 20× dilution is by weight of a
solution whose solids rise with the xylose dose.

### 4c. Functional / structural data (not transcribed; one line each)

- **Table 1 / Table 2** (FTIR, 8 bands × 6 arms, before / after printing): values printed as "Absorption peak
  strength, T (%)", i.e. transmittance (sec. 5, flag 3).
- **Table 3** (XRD 2θ and peak intensity, cps; read on the 150 dpi render): CK has a 5.66° peak (840 cps) that is
  absent in all xylose arms; the authors attribute it to unreacted hydrolysate components [Q]. Structural, not chemical.
- **Table 4** (DSC Tg onset / mid / end, °C): 26.0-28.7 °C in all arms; Fig. 7 heat-flow curves.
- **Fig. 5** photographs of printed hearts (colour by eye: CK greyish-cream, xylose arms cream to beige); **no
  colour measurement**.
- **Fig. 6** cohesion (0.35-0.47) and adhesiveness (~9-24 J) vs xylose dose.

## 5. Flags and internal consistency

1. [!] **No reaction-product chemistry at all.** The paper is titled on a Maillard product but never measures the
   extent of reaction (no A420 / A294 / pH / amine / sugar / volatile). Whether, and how far, the Maillard reaction
   ran in each arm is asserted, not shown [Q].
2. [!] **Dose of xMRP in the dough not stated** (volume of 20×-diluted product per mass of starch); the dose-response
   curves (Figs. 6, 8) are against xylose in the reaction, not against MRP in the dough.
3. [!] **FTIR "absorption peak strength, T (%)"** is transmittance; higher T means weaker absorption, yet the text
   reads rising T as rising intensity ("significantly increased the strength"). (The text's separate statement that
   transmittance fell after printing for 2 g at 3358 / 2932 / 1160 / 991 cm⁻¹ and for 4 / 6 g at 574 cm⁻¹ does agree
   with Tables 1-2, e.g. 2 g at 3358: 46.2 → 38.1.) Does not affect the model.
4. [!] **Fig. 6 caption says "before and after 3D printing"** but each panel shows one curve with no legend; which
   state is plotted is not stated. The text says cohesion increased after printing (so a second series must exist).
5. [!] **Antioxidant sample ambiguity** (sec. 3): methods = "xMRP solution", figure = "materials".
6. [!] **Adhesiveness text vs figure**: text says adhesion "reached the maximum value" at 6 g; Fig. 6 (right) shows
   the **minimum** near 6-8 g (~9 J) and maxima at 2 g (~24 J) and 10 g (~23 J) [F, ±1 J]. The caption itself says
   "first increasing, then decreasing and finally increasing".
7. [!] FTIR peak maximum at **8 g** (Table 1) vs "best" at 6 g everywhere else; the authors note it but do not
   reconcile.
8. [!] Replicate number never stated; no ± in Tables 1-4.
9. [!] Typos: "xMPRs" alongside "xMRPs", "Millard" (abstract), "o-diazophane", "ferrous sulfide"; the text layer
   renders 2θ as "2⊔" (the render shows 2θ correctly).
10. [!] The protein-hydrolysis optimum numbers (40.3 °C, pH 7.02, 6.80 %, 10.6 %) come from an unreported
    optimisation; no reference is given for it.

## 6. Existing records in the repo that cite this paper

`grep -rIn` (2026-11-02) over `data/`, `docs/`, `src/`, `scripts/`, `tasks/` and `results/` for `10327`,
`jsfa.10327`, `zhou2020`, `Quancheng`, `xMRP`, `pea protein enzymatic`, `xylose.pea`, `3D.print`: **no record cites
this paper** (the only `10327` hits are digits inside floating-point numbers in `results/`). `data/keys/papers.yml`
has no row for the DOI. The main checkout's `data/articles/DOWNLOAD_QUEUE.md` (gitignored) does not mention it either. **No
existing dossier cites it**, and no dossier on disk describes a pea-protein-hydrolysate + xylose Maillard pot
(the corpus's hydrolysate + xylose systems are goat viscera, `desousafontes2024`; corn peptide; tilapia head; sesame
meal without added xylose, `ruan2025`; and wheat-protein hydrolysate in `nedvidek1992`). Nothing to check field by
field.

## 7. Proposed ledger entries

**None.** Reasons:
- The model's trunk lane (xylose + amine → Amadori → deoxyosones → furfural / browning) and sulfur lane need at least
  one of: sugar or amine loss, a dicarbonyl / furan / pyrazine / sulfur volatile, A420 / A294, or pH. The paper has
  none (sec. 1, answer 2).
- The amine charge cannot be specified within an order of magnitude (no DH, no free amino N, no soluble-protein
  content), so even a "did it brown?" check is not simulable.
- The quantitative series that exist (Fig. 8 scavenging; FTIR, XRD, TPA, DSC) are properties of a printed starch
  dough with an unknown xMRP loading; no lane or matrix-layer observable corresponds to them.
- Precedent: `huber1990` and `kallen1973` were ingested with "no ledger-ready data" (ledger batch-2 header comments,
  e.g. `trunk.yml` l. 7815).

If the keeper wants the pot design findable anyway, the natural home would be a `trunk.yml` entry with
`kind: other`, `status: text_only`, `quant_class: unknown`, `simulable: false`, `flags: [no_numbers]` (the
`martins2003b_sec324_sugar_controls_not_shown` pattern); I do not recommend it (P(it is ever used) ~2 %).

## 8. Supplementary material

**None referenced.** The text layer has no "Supporting Information", "Appendix", "Table S" or "Fig. S" string; the
article ends with acknowledgements and references. Nothing model-relevant is missing from the PDF on disk; the
missing information (DH, MW distribution, pH, browning) was simply not measured or not reported.
