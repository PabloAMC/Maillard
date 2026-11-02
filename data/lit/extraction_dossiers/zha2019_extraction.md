# Zha, Yang, Rao & Chen 2019 — EXTRACTION (file stem `zha2019`)
DOI: 10.1021/acs.jafc.9b04099. J. Agric. Food Chem. **67**(36):10195-10206 (received 1 Jul 2019, accepted and
published online 22 Aug 2019). Crossref 2026-11-02 (`api.crossref.org/works/10.1021/acs.jafc.9b04099`): title
"Gum Arabic-Mediated Synthesis of Glyco-pea Protein Hydrolysate via Maillard Reaction Improves Solubility, Flavor
Profile, and Functionality of Plant Protein", authors Fengchao Zha, Zhongyu Yang, Jiajia Rao, Bingcan Chen, JAFC 67,
issue 36, pp. 10195-10206, issued 2019-08-22. **Matches the PDF and the brief.**
Stem caveat: **not to be confused with** `zhai2019_extraction.md` (Zhai et al., Jiangnan, DOI 10.1021/acs.jafc.9b04374)
nor with the repo's `zha_2020_ppi_glycation_aggregation` records (Zha et al. 2020, DOI 10.1021/acs.jafc.0c04281, same
first author, not on disk). The group's two earlier conjugation papers, cited here as refs 16 and 18 for every method
detail (Zha, Dong, Rao & Chen 2019 Food Hydrocoll. 92:30, pea protein concentrate + GA; and 2019 Food Chem. 285:130,
pea protein isolate + GA emulsions), are not on disk.

### Commercial pea protein hydrolysate (PURIS 870H, DH 12 %) + gum arabic 1:4 (w/w), hydrated, pH 7.0, freeze-dried, then 5 g incubated in an **open** 60 mL jar in a desiccator at **60 °C / 79 % RH** (saturated KBr) for **0, 1, 3, 5 days**. Gives: free amino groups, A304, A420 and L*/b* at 4 times (with a PPH-alone control), SEC-MALLS Mw / Mn / Rg per fraction (Table 1), six "beany" volatiles by external-standard HS-SPME-GC-MS (Table 2, no PPH-alone incubated control), and hydroperoxide / hexanal time courses in 2 % corn-oil emulsions at 37 °C. **No viscometry, no Mark-Houwink exponent, no sensory panel, no moisture or a_w of the powder, one temperature, one ratio.** The 95-volatile list (pyrazines, thiophenes, Strecker aldehydes) is only in SI Table S1, not on disk.

**Read 2026-11-02 from `data/articles/zha2019.pdf`** (12 pp., ACS, born-digital). Prose from `pdftotext` (plain and
`-layout`) of all pages; all 12 pages rendered at 110 dpi and pp. 5, 6, 7, 10 viewed. **Table 1 (p. 5) read from a
250 dpi crop** (`ingest/zha2019_work/t1hi-05.png`); **Table 2 (p. 6, rotated) from a 220 dpi crop rotated upright**
(`t2rot.png`); the text layer scrambles both. **Fig. 3A/3B/3C (p. 7) and Fig. 6A/6B (p. 10) digitised from 400 dpi crops**
(`f3A-07.png`, `f3Bplot-07.png`, `f3C-07.png`, `f6A-10.png`, `f6B-10.png`), by pixel position against the printed axis
ticks; read errors stated per table. Fig. 1 (SDS-PAGE, FTIR, SEM), Fig. 2 (SEC traces), Fig. 4 (PCA, chromatograms) and
Fig. 5 (solubility, emulsions) viewed at 110 dpi only, not digitised (qualitative or functional).
All scratch files are under
`/private/tmp/claude-501/-Users-pabloantoniomorenocasares-Developer-Maillard/7ae75f75-b2b3-4540-b5f1-5f645399920f/scratchpad/ingest/zha2019_work/`.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Plant Sciences (Zha, Rao, Chen) and Department of Chemistry and Biochemistry (Yang),
**North Dakota State University, Fargo**. Corresponding authors **Bingcan Chen** (Plant Sciences; senior, last author)
and **Zhongyu Yang** (Chemistry; the SEC-MALLS side). **Not in `data/keys/laboratories.yml`: NEW.** Proposed id
**`ndsu_chen_rao`**, string "North Dakota State Univ. (Chen B / Rao J)". Evidence for one group: Chen and Rao co-author
this paper and `gao2020_extraction.md` (Gao, Shen, Lan, Cui, Ohm, Chen & Rao 2020 Food Res Int 131:109045; Rao
corresponding), same department, same six-marker "beany" HS-SPME panel (hexanal, 1-octen-3-ol, 1-octen-3-one,
2-pentylfuran, 3-methyl-1-butanol; here acetophenone instead of 1-pentanol). The three gao2020 ledger `lab` strings
("North Dakota State University / USDA-ARS Fargo (Rao)" in lipid.yml, "North Dakota State Univ. / USDA-ARS Fargo (Rao)"
in matrix.yml, "North Dakota State (Rao J.)" in aroma_other.yml) are currently unmapped and would belong to the same
group (confidence ~80 % that Chen and Rao share methods and instruments; the conservative leave-one-lab-out choice is
to merge). Same-lab dossier: `gao2020_extraction.md` only (grep of `extraction_dossiers/` for North Dakota / NDSU /
Bingcan, 2026-11-02).

## 1. What this paper answers for the model

| model need | what Zha 2019 gives | tag |
|---|---|---|
| amine loss vs time in a pea-protein + polysaccharide system | Fig. 3B: free amino groups (mmol/g protein) at 0/1/3/5 d, PPH alone and PPH+GA, 60 °C / 79 % RH; text: −26.5 % (PPH+GA) and −9.8 % (PPH) at 3 d | [F] / [M] |
| assay behind the amine numbers | **not stated** (spectrophotometer only; "following our previous studies"); OPA vs TNBS unknown | [!] |
| Amadori / intermediate (A294-type) | **A304** (not A294), Fig. 3A, same 4 times, both arms; sample concentration not stated | [F] |
| browning | A420 (Fig. 3A) and L*, b* (Fig. 3C; a* only as bubble size), same 4 times | [F] |
| temperature series / Ea | **none**: 60 °C only | — |
| water activity / moisture | RH 79 % set by saturated KBr; **powder moisture and equilibration not measured**; jar open | [M] / [!] |
| sugar inventory (reducing ends, free sugars) | **none**; GA "~90 %" GA; PPH carbohydrate ~6 % (supplier spec) | [M] |
| conjugate size | Table 1 SEC-MALLS Mw, Mn, Rg, Mw/Mn, mass fraction per fraction; N ≈ 1.2 mol GA / mol PPH | [M] |
| beany (lipid) volatiles before / after | Table 2: hexanal, 1-octen-3-ol, 1-octen-3-one, 2-pentylfuran (+ 3-methyl-1-butanol, acetophenone), external standard, 0/1/3/5 d; **units contradictory** (sec. 5) | [M] |
| Maillard volatiles (pyrazines, furans, thiophenes, Strecker aldehydes) | named in text; numbers only in SI Table S1 (relative to E-2-octenal), **not on disk** | [Q] |
| control for volatilisation (heated PPH alone, closed jar) | **absent** for volatiles | [!] |
| hexanal formation in emulsions | Fig. 6B: hexanal (ng/mL) 0-10 d at 37 °C in 2 % corn-oil emulsions with PPH / mixture / glyco-PPH-1; Fig. 6A LOOH (µM) | [F] |
| sensory | **none** (PCA of volatiles and a cited pea sensory study only) | [C] |
| Mark-Houwink exponent, viscosity | **none**: no viscometer on the SEC line (UV, RI, MALLS only) | [!] |

**Answers.**

1. **Maillard conditions and time course** (confidence ~95 % on the conditions as printed). PPH : GA **1 : 4 by mass**,
   hydrated 1 : 2 (w/v) in deionised water 24 h at 22 °C, pH adjusted to 7.0, freeze-dried; **dry-state**: 5 g of the
   lyophilised mix in an uncovered 60 mL glass jar on a perforated plate in a desiccator held at **79 % RH (saturated KBr)
   inside a 60 °C incubator**, for **0, 1, 3, 5 days** (Methods p. 2). PPH alone was incubated the same way for SDS-PAGE,
   A304/A420, free amino groups, colour and solubility. **Time courses exist** for free amino groups (Fig. 3B), A304 and
   A420 (Fig. 3A) and colour (Fig. 3C), 4 points each, both arms; **no degree-of-graft measure** other than the
   SEC-MALLS N and the free-amino loss, and the amine assay is not named (refs 16, 18). For the trunk lane these are the
   only numbers, and **they are weak** (P(usable as a fit table) ~10 %; P(usable as a directional check: GA roughly
   triples the amine loss of heated PPH, most of it in day 1, then a plateau) ~60 %). Reasons: (a) the free-amino level,
   **~6.2 mmol/g protein at t = 0**, is ~4 × what a 12 %-hydrolysed pea protein can carry (lysine ~0.5 mmol/g plus
   ~1 mmol/g new α-amino at DH 12 %) [Z, ~80 % confidence of a unit or calibration problem]; (b) **the extra amine loss
   due to GA (~1.06 mmol/g protein at 3 d) is ~180 × the reducing ends GA can supply** (GA number-average M from
   Table 1 ≈ 7.8 × 10^5 g/mol → ~1.3 µmol reducing ends per g GA → ~6 µmol per g protein in a 1 : 4 mix) [Z]; so the
   amine loss cannot be GA reducing-end glycation, and must come from free sugars / GA hydrolysis products, lipid
   carbonyls or assay interference (~85 % confidence); (c) a_w is a set point, not a measured state; (d) one T.
2. **Volatiles** (confidence ~90 % on what is printed, ~60 % on the interpretation below). Table 2 quantifies six
   markers by external-standard calibration (LOD 0.008-0.025 mg/L, recovery 91.8-98.0 %): hexanal, 3-methyl-1-butanol,
   2-pentylfuran, 1-octen-3-one, 1-octen-3-ol, acetophenone, in PPH, GA, the 0-d mixture and glyco-PPH 1/3/5 d. **Units
   contradict**: footnote "milligrams per gram of dry weight" vs text "hexanal (3.54 ppm)"; ppm (µg/g) is the plausible
   one (~90 %). Hexanal falls 3.17 → 1.28 → 0.78 → 0.42 (mixture 0 d → 1 → 3 → 5 d). Maillard volatiles (7 pyrazines incl.
   2-ethyl-3,5-dimethyl- and 3-ethyl-2,5-dimethylpyrazine; 2-methyl- and 2-ethylthiophene; benzaldehyde, furfural,
   5-methylfurfural; 4,5-dihydro-5-propyl/butyl-2(3H)-furanone at 5 d) are named [Q] but quantified only in SI Table S1
   (relative to E-2-octenal), not on disk. **Attribution**: the authors attribute the beany reduction to conformational
   change of PPH "resulting in the release of hydrophobic beany flavors" [V, p. 9] — i.e. **loss by release, not masking
   or binding and not consumption** [P]; they also call glyco-PPH aroma "nonspecific" from Maillard products [Q].
   **The data cannot tell the mechanisms apart**: no heated-PPH-alone volatile arm, an open jar at 60 °C for days, no
   headspace / mass balance. One internal clue [Z]: **all six markers drop by a similar fraction in day 1** (glyco-1 /
   mixture-0 = 0.31-0.40, including the alcohols 3-methyl-1-butanol and 1-octen-3-ol and the ketone acetophenone, which
   cannot react with amines), which favours a physical loss (stripping from the open jar, or a matrix change that lowers
   SPME recovery) over carbonyl-amine consumption. My split: P(volatilisation / stripping dominant) ~55 %, P(changed
   matrix retention or extraction efficiency dominant) ~30 %, P(Maillard consumption of hexanal / octenone dominant)
   ~15 %. Formation of lipid volatiles during incubation (PPH has ~8 % lipid) runs the other way and is not separable.
3. **Sensory**: none. Only an OAV argument with cited thresholds "40−164 ppb" [C] and a PCA of volatiles.
4. **Proposed ledger entries**: seven (sec. 7): trunk (free amino groups; A304 / A420; colour), aroma_other (SEC-MALLS
   molar mass; Table 2 non-lipid rows; SI Table S1, supplement_missing), lipid (Table 2 lipid rows; Fig. 6 emulsion LOOH
   / hexanal). None is simulable today.
5. **The repo's three records** (sec. 6): **none of their eight numbers is in the paper** (confidence ~97 %). They trace
   to a Gemini deep-research report whose reference **21** is this DOI; the "Ref. 21" in the record ids is that report's
   footnote number (in the paper itself ref. 21 is Bradford 1976).

## 2. Conditions

- **Proteins** [M]: PURIS pea protein 870 (PPI) and **PURIS 870H (PPH, its hydrolysate)**, Cargill. PPH spec: protein
  ~80 %, moisture ~6 %, ash ~5 %, carbohydrate ~6 %, lipid ~8 %. Degree of hydrolysis 12 % (p. 6 text). Enzyme not
  stated. Used as received.
- **Polysaccharide** [M]: TIC Pretested gum arabic spray-dried powder (TIC Gums): moisture ~6.7 %, GA ~90 %, protein
  ~3.2 %, minerals ~0.2 %. Free-sugar content not given.
- **Mixture** [M]: PPH : GA **1 : 4 (w/w)**; hydrated in deionised water 1 : 2 (w/v) 24 h, 300 rpm, 22 °C; **pH adjusted
  to 7.0** (reagent not stated); lyophilised. pH after incubation **not measured**.
- **Heating** [M]: 5 g powder in a VWR 60 mL straight-sided jar, **uncovered**, on a perforated porcelain plate in a
  desiccator over saturated KBr (**79 % RH**) inside a preheated incubator at **60 °C** (Heratherm IMH180). Times **0, 1,
  3, 5 days**. Equilibration time and moisture of the powder not reported [!]. Headspace = air; O2 not controlled.
- **Controls** [M]: PPH alone through the same incubation (SDS-PAGE, A304, A420, amines, colour, solubility); the 0-d
  mixture "PPH-GA-0"; GA alone (volatiles, SEC).
- **Emulsions** (Fig. 5-6) [M]: 2 wt % Mazola corn oil (SFA ~14.3 %, MUFA ~28.6 %, PUFA ~57.1 %) in 98 wt % of 0.20 wt %
  emulsifier in 10 mM PBS pH 7.0; blender 2 min then APV-Gaulin two-stage homogeniser 5000 / 500 psi × 3; 0.04 % sodium
  azide; storage at **37 °C**, 0-10 d. Only glyco-PPH-1 (1 d conjugate) tested against PPH and the 0-d mixture.
- **Replicates** [M]: "at least two independent experiments" for sample preparation; measurements in triplicate;
  ± = SD of triplicates; one-way ANOVA, Tukey p < 0.05 (SAS 9.4). Table 2 states n = 3.

## 3. Analytical methods

- **Free amino groups** [M/!]: "recorded using a Shimadzu UV-1100 spectrophotometer"; reagent, wavelength, standard
  and sample concentration **not given** (deferred to refs 16, 18). Units mmol/g protein (Fig. 3B axis).
- **A304 / A420** [M]: absorbance of glyco-PPH "solution" at 304 nm (taken as Amadori compounds) and 420 nm
  (melanoidins), refs 16, 33; concentration and path not given. **Note: 304 nm, not 294 nm.**
- **Colour** [M]: Minolta CR-310 chroma meter, L*, a*, b*.
- **SEC-MALLS** [M]: 0.10 g sample in 10 mL 10 mM PBS pH 7.0, 2 h hydration, 2000 g 30 min, 0.45 µm nylon filter (the
  insoluble part is discarded before analysis [!]); Agilent 1200, polySep-GFC-P guard + linear column, 0.3 mL/min, 100 µL;
  DAD 280 nm, RI (Agilent 1362A), Wyatt DAWN HELEOS II (λ 661 nm [sic; HELEOS II lasers are 658 nm]); **one dn/dc =
  0.141 mL/g ("gum arabic") applied to every sample, including PPI and PPH** [!]; RI calibrated with NaCl; ASTRA 7.1.2.5;
  22 °C. **No viscometer** → no intrinsic viscosity, no Mark-Houwink exponent. N (GA per PPH) by Eq. 1-2:
  N = (M̄w2 − M̄w1)/M̄w3 with M̄w = Σ(Fi Mi)/Σ Fi (PPH, glyco-PPH, GA).
- **Volatiles (powders)** [M]: HS-SPME (CTC), Agilent 7890B / 5977A, NIST identification; parameters "in our previous
  work" (ref. 16). Table 2: six markers with authentic standards ("St"), quantitative / qualifier ions, external
  calibration with LOD (mg/L; the footnote calls it "limit of quantitation"), recovery %, RSD %. PCA and Table S1 use
  relative content normalised to **E-2-octenal as internal standard**. Sample mass per vial not stated here.
- **Emulsion oxidation** [M]: lipid hydroperoxides by the Chen-McClements-Decker method (ref. 22); hexanal by the
  method of "Zhao and co-workers" [sic; ref. 18 is Zha et al.], calibration with an authentic standard, **LOD 7.89 ng/mL**.
  Lag phase = time of a sudden increase.
- Solubility (Bradford), particle size (Mastersizer 3000, d43), ζ-potential (Nano-ZS): functional, not extracted.

## 4. Data

### 4a. Table 1 (p. 5): SEC-MALLS molecular parameters [M] (± as printed, relative %)

| sample | fr. | mass fr. (%) | Mw (g/mol) | Mn (g/mol) | Rg (nm) | I = Mw/Mn | mass recovery (%) |
|---|---|---|---|---|---|---|---|
| PPI | 1 | 48.7 | 3.59 × 10^5 ± 12.6 % | 1.29 × 10^5 ± 6.41 % | 34.0 ± 0.81 % | 2.76 ± 14.2 % | 98.2 |
| PPI | 2 | 44.6 | 1.72 × 10^5 ± 8.42 % | 5.95 × 10^4 ± 8.60 % | 28.6 ± 1.82 % | 2.89 ± 2.14 % | 11 |
| PPI | 3 | 6.70 | 6.01 × 10^4 ± 3.41 % | 2.54 × 10^4 ± 4.11 % | 22.4 ± 0.90 % | 2.36 ± 3.07 % | 7.1 |
| PPH | 1 | 58.9 | 2.81 × 10^5 ± 5.20 % | 8.19 × 10^4 ± 6.60 % | 32.8 ± 0.10 % | 3.43 ± 8.37 % | 94.1 |
| PPH | 2 | 34.2 | 1.52 × 10^5 ± 0.82 % | 4.16 × 10^4 ± 8.60 % | 22.6 ± 1.80 % | 3.69 ± 28.8 % | 34 |
| PPH | 3 | 6.90 | 5.61 × 10^4 ± 1.11 % | 1.98 × 10^4 ± 4.12 % | 20.4 ± 1.93 % | 2.83 ± 7.36 % | 12.1 |
| GA | 1 | 80.1 | 2.17 × 10^6 ± 0.48 % | 1.72 × 10^6 ± 0.43 % | 30.5 ± 0.15 % | 1.26 ± 0.64 % | 86.3 |
| GA | 2 | 19.9 | 2.64 × 10^5 ± 5.49 % | 2.43 × 10^5 ± 3.84 % | 27.7 ± 0.72 % | 1.08 ± 6.70 % | 19.3 |
| PPH-GA-0 | 1 | 31.5 | 2.48 × 10^6 ± 0.39 % | 2.41 × 10^6 ± 0.38 % | 36.8 ± 1.19 % | 1.03 ± 0.54 % | 92.5 |
| PPH-GA-0 | 2 | 68.5 | 3.23 × 10^5 ± 1.34 % | 5.18 × 10^4 ± 0.80 % | 26.3 ± 1.52 % | 6.34 ± 1.56 % | 89.5 |
| glyco-PPH-1 | 1 | 47.4 | 4.17 × 10^6 ± 0.38 % | 2.87 × 10^6 ± 0.36 % | 38.8 ± 0.18 % | 1.45 ± 0.52 % | 99.1 |
| glyco-PPH-1 | 2 | 52.6 | 6.22 × 10^5 ± 5.37 % | 1.31 × 10^5 ± 1.39 % | 34.2 ± 0.21 % | 4.75 ± 5.55 % | 90.1 |
| glyco-PPH-3 | 1 | 39.9 | 3.66 × 10^6 ± 1.19 % | 3.27 × 10^6 ± 0.39 % | 32.2 ± 0.18 % | 1.14 ± 0.61 % | 99.4 |
| glyco-PPH-3 | 2 | 60.1 | 7.75 × 10^5 ± 7.29 % | 3.97 × 10^5 ± 6.67 % | 28.5 ± 0.33 % | 1.95 ± 9.88 % | 68.7 |
| glyco-PPH-5 | 1 | 36.9 | 3.36 × 10^6 ± 0.38 % | 3.33 × 10^6 ± 0.35 % | 36.7 ± 0.14 % | 1.01 ± 0.52 % | 98.1 |
| glyco-PPH-5 | 2 | 63.1 | 7.96 × 10^5 ± 6.64 % | 4.68 × 10^5 ± 5.33 % | 26.8 ± 0.25 % | 1.70 ± 8.52 % | 78.6 |

Text [M]: PPI fractions ~340 kDa (11S hexamer), ~160 kDa (7S trimer), ~15 kDa (2S dimer) [C]; PPH fractions ~254,
~158, ~13 kDa and "MW ≈ 280 kDa" for a legumin-like fraction; GA AG ~80 % at ~2.17 × 10^6, AGP ~20 % at ~2.64 × 10^5;
GA Rg "30.5 nm for the AG components and 27.7 nm" for AGP; **N ≈ 1.2 mol GA per mol PPH after 1 d**; fraction 1 share
47.4 → 40 → 37 % at 1 → 3 → 5 d.

Derived [Z]:
- Mass-weighted Mw (Eq. 2): PPH 2.21 × 10^5, GA 1.79 × 10^6, PPH-GA-0 1.00 × 10^6, glyco-PPH-1 2.30 × 10^6,
  glyco-PPH-3 1.93 × 10^6, glyco-PPH-5 1.74 × 10^6 g/mol. **N = (2.30 − 0.221)/1.79 = 1.16** for glyco-PPH-1 (reproduces
  the printed ~1.2). **The same formula applied to the unheated mixture gives N = 0.44**: the method does not subtract
  the free GA (80 % of the sample mass), so N is not a graft ratio.
- Stoichiometric ceiling: number-average M (harmonic over fractions) PPH 5.3 × 10^4, GA 7.8 × 10^5 g/mol → in a 1 : 4
  mix there are ~0.27 mol GA per mol PPH chains. **N = 1.2 exceeds the GA available** (~0.27 on a number basis; ~0.5
  on a weight-average basis) [!].
- Direction of the conjugate: fraction 1 Mw 2.48 → 4.17 MDa at 1 d (+68 %), then 3.66, 3.36 MDa; fraction 2 Mw 0.32 →
  0.62 → 0.78 → 0.80 MDa. Rg fraction 1: 36.8 → 38.8 → 32.2 → 36.7 nm (non-monotonic). Mw/Mn fraction 1: 1.03 → 1.45 →
  1.14 → 1.01; fraction 2: 6.34 → 4.75 → 1.95 → 1.70 (**narrows** with time).
- Mw/Mn recomputed from Mw and Mn: agrees with the I column to ±0.02 except PPH-GA-0 fraction 2 (6.24 vs 6.34).

### 4b. Table 2 (p. 6): selected beany volatiles in the powders, mean ± SD (n = 3) [M]
Units as printed: "milligrams per gram of dry weight" (footnote a); text gives the PPH column as **ppm** (sec. 5, flag 1).
Letters = Tukey groups across samples. GA column: nd for all six.

| no. | compound | quant. ion | PPH | PPH-GA-0 | glyco-PPH-1 | glyco-PPH-3 | glyco-PPH-5 | LOD (mg/L) | recovery % | RSD % |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | hexanal | 56 | 3.54 ± 0.04 a | 3.17 ± 0.06 a | 1.28 ± 0.07 b | 0.78 ± 0.02 c | 0.42 ± 0.04 d | 0.025 | 92.59 | 4.68 |
| 2 | 3-methyl-1-butanol | 55 | 0.45 ± 0.02 a | 0.36 ± 0.01 b | 0.12 ± 0.01 c | 0.05 ± 0.01 c | 0.03 ± 0.01 c | 0.023 | 96.40 | 1.64 |
| 3 | 2-pentylfuran | 81 | 0.19 ± 0.04 a | 0.16 ± 0.02 a | 0.05 ± 0.01 b | 0.02 ± 0.01 c | 0.01 ± 0.00 c | 0.008 | 91.83 | 2.61 |
| 4 | 1-octen-3-one | 55 | 1.07 ± 0.01 a | 0.82 ± 0.02 a | 0.33 ± 0.01 b | 0.09 ± 0.01 c | 0.04 ± 0.02 d | 0.015 | 97.97 | 7.12 |
| 5 | 1-octen-3-ol | 57 | 1.10 ± 0.03 a | 0.95 ± 0.02 a | 0.34 ± 0.02 b | 0.19 ± 0.01 c | 0.11 ± 0.01 d | 0.016 | 95.86 | 0.12 |
| 6 | acetophenone | 105 | 0.40 ± 0.07 a | 0.45 ± 0.03 a | 0.17 ± 0.04 b | 0.13 ± 0.01 b | 0.07 ± 0.01 c | 0.012 | 94.83 | 6.14 |

Qualifier ions [M]: 56, 44, 41 / 55, 70, 42 / 81, 138, 53 / 55, 70, 27 / 57, 72, 43 / 105, 120, 77. The PPH column is
**unincubated PPH** (no heated-PPH-alone column exists).

Derived fractions remaining vs PPH-GA-0 [Z] (1 / 3 / 5 d): hexanal 0.40 / 0.25 / 0.13; 3-methyl-1-butanol 0.33 / 0.14 /
0.08; 2-pentylfuran 0.31 / 0.13 / 0.06; 1-octen-3-one 0.40 / 0.11 / 0.05; 1-octen-3-ol 0.36 / 0.20 / 0.12;
acetophenone 0.38 / 0.29 / 0.16. Apparent first-order hexanal loss: 0.91 d⁻¹ over day 0-1, 0.25 d⁻¹ over 1-3, 0.31 d⁻¹
over 3-5 (biphasic).

### 4c. Fig. 3B: free amino groups (mmol/g protein), 60 °C / 79 % RH [F] (read ±0.03; error bars ±0.10-0.26)

| day | PPH alone | PPH + GA (glyco-PPH) |
|---|---|---|
| 0 | 6.16 | 6.28 |
| 1 | 5.66 | 4.84 |
| 3 | 5.55 | 4.61 |
| 5 | 5.43 | 4.43 |

Text [M]: −26.5 % at 3 d (PPH+GA), −9.8 % at 3 d (PPH alone). Digitised losses [Z]: PPH+GA 22.9 / 26.6 / 29.5 % at
1 / 3 / 5 d (reproduces 26.5); PPH 8.1 / 9.9 / 11.9 % (reproduces 9.8). Differential (GA-attributable) loss at 3 d
≈ 1.06 mmol/g protein [Z]. Curves in the figure are hand-drawn smooth lines, not fits.

### 4d. Fig. 3A: A304 and A420 [F] (A304 read ±0.008, A420 read ±0.004; two y-axes: A304 0-0.40, A420 0-0.24)

| day | A304 PPH | A304 glyco-PPH | A420 PPH | A420 glyco-PPH |
|---|---|---|---|---|
| 0 | 0.051 | 0.016 | 0.007 | 0.004 |
| 1 | 0.113 | 0.258 | 0.033 | 0.110 |
| 3 | 0.120 | 0.288 | 0.036 | 0.148 |
| 5 | 0.125 | 0.211 | 0.045 | 0.182 |

A304 (glyco) peaks at 3 d and falls by 5 d; A420 keeps rising [F/Q]. Error bars on glyco-PPH ±0.01-0.03.

### 4e. Fig. 3C: colour [F] (L* read ±0.4, b* read ±0.5; a* encoded only as bubble size, not readable)

| day | PPH L* | PPH b* | glyco-PPH L* | glyco-PPH b* |
|---|---|---|---|---|
| 0 | 82.8 | 22.0 | 83.6 | ~14.3 |
| 1 | 79.0 | 22.5 | 74.5 | 24.4 |
| 3 | 78.5 | 23.5 | 72.6 | 27.6 |
| 5 | 78.1 | 23.9 | 72.4 | 28.7 |

### 4f. Fig. 6: 2 % corn-oil-in-water emulsions, pH 7.0, 37 °C [F] (x axis mislabelled "Conjugation time/day"; it is storage time)

LOOH (µM, read ±0.4); 0 at days 0-3 for all three except as noted:

| day | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| PPH | 1.7 | 3.0 | 7.4 | 9.2 | 11.0 | 16.5 | 21.5 |
| PPH-GA-0 | 1.4 | 2.5 | 4.8 | 7.6 | 9.0 | 15.5 | 19.5 |
| glyco-PPH-1 | ~0 | ~0.1 | 0.3 | 2.2 | 5.8 | 10.8 | 17.4 |

Hexanal (ng/mL, read ±3; LOD 7.89 ng/mL):

| day | 0-3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| PPH | 0 | ~3 | 17 | 17 | 41 | 41 | 107 | 127 |
| PPH-GA-0 | 0 | 0 | 0 | 0 | ~2 | 14 | 33 | 39 |
| glyco-PPH-1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 13 |

Text [Q]: glyco-PPH hexanal "lower than LOD even after 9 days" (consistent: day 10 ≈ 13 ± 7 ng/mL is the first point above
LOD). Lag phases read [Z]: hexanal ~4 d (PPH), ~7 d (mixture), ~9 d (glyco-PPH-1).

### 4g. Functional data (not extracted)
Solubility at pH 7.0 19.4 → 26.2 % after 1 d conjugation, then falls (Fig. 5A); emulsion d43 20.7 (PPH), 5.08 (mixture),
0.75 µm (glyco-PPH-1) (Fig. 5B); IEP ~3.5 → 2.5; ζ −36.5 mV (PPH, pH 6) (Fig. 5C) [M]. SDS-PAGE: new band >250 kDa and a
~37 kDa "PPH dimer"; FTIR new bands 1592, 1673, 3335 cm⁻¹; SEM morphology (Fig. 1) [Q].

## 5. Flags and internal consistency

1. [!] **Table 2 units**: footnote "mg per g of dry weight" vs text "3.54 ppm" for the same hexanal cell. mg/g would be
   3540 mg/kg hexanal in a protein powder (implausible; the corpus's pea-isolate hexanal is 15-83 mg/kg); µg/g (ppm) is
   the plausible reading (~90 %). LOD is in mg/L of an unstated solution; the footnote calls LOD a "limit of
   quantitation".
2. [!] **Table 2 basis**: PPH-GA-0 is 20 % PPH by mass, yet its hexanal (3.17) is 90 % of pure PPH (3.54). If per g of
   sample dry weight, the 1 : 4 dilution predicts ~0.71 [Z]. So either the basis is per g PPH / equal protein (not stated;
   ~65 %), or freeze-drying the mixture did not dilute the signal for some other reason. The comparison the authors make
   (glyco-PPH vs PPH, "more than 2-fold less") mixes this basis question with the conjugation effect; glyco vs PPH-GA-0
   is the fair contrast (sec. 4b).
3. [!] **No heated-PPH-alone volatile arm and an uncovered jar** at 60 °C for up to 5 d: volatilisation, oxidation and
   Maillard consumption are confounded (sec. 1, ans. 2).
4. [!] **Free amino groups ~6.2 mmol/g protein** at t = 0, ~4 × the plausible ceiling for DH 12 % pea protein [Z]; assay
   not named.
5. [!] **Amine loss vs GA reducing ends**: the GA-attributable loss (~1.06 mmol/g protein) is ~180 × the GA reducing ends
   in the mix [Z]; the authors still read it as "cross-linking of GA to PPH" [Q].
6. [!] **N ≈ 1.2 mol GA/mol PPH** is above the ~0.27 mol GA per mol PPH chains present, and the same formula gives 0.44 for
   the unheated mixture [Z]. Not a graft ratio.
7. [!] **One dn/dc (0.141 mL/g, gum arabic) for all samples**, including the pure proteins (protein dn/dc is ~0.185);
   PPI / PPH Mw would be biased ~30 % high [Z, ~70 % confidence in the direction]. Insoluble material removed by
   centrifugation and filtration before SEC.
8. [!] **"Mass recovery (%)" per fraction sums to > 100 %** within a sample (PPI 116, PPH 140, PPH-GA-0 182, glyco-PPH-1
   189 %) [Z]; the column's meaning is undefined.
9. [!] Text cites "(Table 2)" for a molar-mass increase (should be Table 1) and mixes units ("16.6 min" vs "17.2 mL").
10. [!] SEM text says "after 6 days conjugation (Figure 1C-6)"; the caption and design have 5 days.
11. [!] Fig. 6 x axis labelled "Conjugation time/day"; it is emulsion storage time at 37 °C.
12. [!] The paper's "Amadori" absorbance is **A304**; A304 and A420 sample concentration not stated, so absolute values
    are not transferable.
13. [!] Hexanal emulsion method credited to "Zhao and co-workers" with ref. 18, which is Zha et al. (typo).
14. [!] MALLS laser stated as 661 nm (HELEOS II uses 658 nm; trivial).
15. [!] a_w is a set point (79 % RH), not a measured powder state; no moisture content after incubation.

## 6. Existing records in the repo that cite this paper

Grep 2026-11-02 (`grep -rIn` of `data/ docs/ src/ scripts/ tasks/` and `results/` for `9b04099`, `jafc_2019_ref21`,
`zha2019`, "gum arabic"): three data records, plus `data/keys/papers.yml` (paper_id `jafc_2019_ref21_pea_gum_arabic_
architecture_anchor`, `dossier: null`, `verification_status: autogenerated_unverified`), the generated
`results/literature/slr_family_reports/slr_family_16_melanoidin_polymerization.md` (scores it 2/8, "Rejected", yet
status `encoded`), `results/validation/key_value_trace_report.md` l. 88-89, `results/validation/citation_verification_
ledger.md` l. 572 ("MATCH ... as cited" — DOI identity only), and the source text
`data/research_corpus/raw/Gemini_deep_research_F.md` (and `_11.md` l. 169, reference list only).

**Where "Ref. 21" comes from.** `Gemini_deep_research_F.md` ("Kinetic Formation, Structural Architecture, and Physical
Impacts of High-Molecular-Weight Melanoidins in Extruded Plant Proteins") footnote **21** (l. 119) is this paper's ACS
URL. Its l. 22-23 and table rows l. 31-32 attribute to "21" the sentence that unconjugated protein fractions of ~290 kDa
and PDI 1.35 become 340-359 kDa, PDI 1.74, Rg 29.5 → 30.5-33.5 nm, and Mark-Houwink 'a' 0.53 → 0.40, described as
"carbohydrate byproducts" conjugation. The repo copied these numbers. **In the paper itself ref. 21 is Bradford 1976**
(protein assay). None of the eight numbers appears in the paper's text, Table 1, figures or SI caption (Table S1 is a
volatile list). Near-coincidences only: PPI (not conjugate) fraction 1 Mw = 3.59 × 10^5 ("359"), PPI hexamer "~340 kDa",
GA AG Rg 30.5 nm, glyco-PPH-5 fraction 2 Mw/Mn 1.70. Confidence ~97 % that the record numbers are not from this paper
(the 3 % covers a figure or SI detail I could not see; the SI is not on disk, but its caption is volatiles only).

Code readers: no module in `src/` or `scripts/` reads the fields `native_mw_kda`, `conjugated_*`, `*_mark_houwink_
exponent`, `*_pdi`, the section `melanoidin_trapping_priors`, or the kind `melanoidin_architecture_state` (grep
2026-11-02). `src/matrix_prior_registry.py` loads the whole priors file; `src/deep_research_runtime_queue.py`,
`src/extrusion_external_closure.py` and `src/extrusion_benchmark_landing.py` load `process_state_calibrations.json`
by id for other records. So **no shipped number depends on these records** (~95 %); they feed reports only.

### 6a. `data/lit/process_state_calibrations.json` → `jafc_2019_ref21_pea_gum_arabic_architecture_state`

| field | repo value | paper | verdict |
|---|---|---|---|
| kind | melanoidin_architecture_state | SEC-MALLS of PPH-GA glycoconjugate mixtures; "melanoidins" appear only as the A420 reading | **mislabelled**: no melanoidin was isolated or sized |
| protein_type | pea_iso | pea protein **hydrolysate** (PURIS 870H, DH 12 %) | **mismatch** |
| source_citation / doi | Zha et al. (2019) / 10.1021/acs.jafc.9b04099 | same | **match** |
| provenance_tier | direct_measurement | the numbers are not in the paper | **not supported** |
| conditions.matrix_label | pea_protein_hydrolysate_gum_arabic_conjugate | PPH : GA 1 : 4 | **match** (ratio missing) |
| conditions.process_state | heated_matrix | dry-state incubation, 60 °C, 79 % RH, 1-5 d, open jar | **imprecise**; T, RH, time, ratio absent |
| conditions.measurement_mode | hpsec_mals_mark_houwink | SEC-MALLS with UV + RI; **no viscometer** | **"mark_houwink" not in paper** |
| native_mw_kda | 290.0 | PPH fractions 281 / 152 / 56 kDa (mass-weighted 221 kDa); text "≈ 280 kDa"; mixture PPH-GA-0 2.48 MDa / 323 kDa | **not in paper** |
| conjugated_mw_kda_range | [340, 359] | glyco-PPH fraction 1: 3360-4170 kDa; fraction 2: 622-796 kDa (1-5 d) | **mismatch, ~10 × too low**; 359 = PPI (unhydrolysed, unconjugated) fraction 1 |
| native_radius_of_gyration_nm | 29.5 | PPH 32.8 / 22.6 / 20.4; PPH-GA-0 36.8 / 26.3 | **not in paper** |
| conjugated_radius_of_gyration_nm_range | [30.5, 33.5] | glyco-PPH fr. 1: 38.8 / 32.2 / 36.7; fr. 2: 34.2 / 28.5 / 26.8 (1 / 3 / 5 d); 30.5 is GA's AG fraction | **not in paper** |
| native_mark_houwink_exponent | 0.53 | not measured | **not in paper** |
| conjugated_mark_houwink_exponent | 0.40 | not measured | **not in paper** |
| native_pdi | 1.35 | PPH Mw/Mn 3.43 / 3.69 / 2.83; PPH-GA-0 1.03 / 6.34 | **not in paper** |
| conjugated_pdi | 1.74 | glyco-PPH 1.45 / 4.75 (1 d), 1.14 / 1.95 (3 d), 1.01 / 1.70 (5 d) | **not in paper** (nearest 1.70); and PDI *narrows* with time, opposite to the "broadening" story |
| what_it_supports [0] ("polymer-architecture shift") | — | an Mw rise at 1 d then a fall (fraction 1 4.17 → 3.36 MDa) is shown, in a mixture dominated by free GA | **partly**: direction at 1 d only |
| what_it_supports [1] ("polymer growth can increase trapping pressure") | — | no binding, retention or headspace partition measured; the authors propose conjugation **releases** hydrophobic beany volatiles [P] | **absent; the paper's own hypothesis runs the other way** |
| comment ("compaction cues") | — | no compaction measure (no Mark-Houwink, no Rg/Mw scaling analysis) | **not in paper** |

### 6b. `data/lit/computational_priors.json` → `melanoidin_trapping_priors` → `jafc_2019_ref21_pea_gum_arabic_architecture_v1`

Same eight numbers and the same reference_conditions as 6a: **every number not in paper; `hpsec_mals_mark_houwink`
not in paper** (verdicts as in 6a). Further fields:

| field | repo value | paper | verdict |
|---|---|---|---|
| source | Zha et al. (2019) | Zha, Yang, Rao & Chen 2019 | **match** |
| provenance_tier / confidence_tier | literature_derived_transfer / medium_low | — | tier labels do not say the numbers are unsourced; should be withdrawn, not tiered |
| section name | melanoidin_trapping_priors | no trapping measured | **not in paper** |
| notes ("trapping-risk direction") | — | no trapping data; release proposed | **absent / opposite** |
| doi_propagation.crossref_verified | title, authors, year 2019, JAFC | Crossref 2026-11-02 agrees | **match** (identity only) |

### 6c. `data/lit/benchmark_intake_registry.json` → `jafc_2019_ref21_pea_gum_arabic_architecture_anchor`

| field | repo value | paper | verdict |
|---|---|---|---|
| citation | "J. Agric. Food Chem. 2019 (Ref. 21)" | Zha et al. 2019 JAFC 67:10195 | **incomplete**; "Ref. 21" is the Gemini report's footnote, not a property of the paper |
| doi | 10.1021/acs.jafc.9b04099 | same | **match** |
| chemistry_family / slr_family_source | melanoidin_polymerization / 16 | protein-polysaccharide glycoconjugation (early Maillard), A420 browning as a side marker | **mislabelled** |
| observable_panel_tags | process_state_calibration, melanoidin_burden, polymer_architecture, hydrolysate | hydrolysate yes; melanoidin burden not measured (A420 only) | **partly** |
| supporting_families | 05, 06 | — | repo decision |
| key_values (eight numbers) | as 6a | — | **all not in paper** |
| what_it_supports [1] ("shift ... toward higher trapping pressure") | — | — | **absent** |
| eligibility / status | runtime_eligible / encoded | — | not justified (no number traceable) |
| requires_primary_data | false | — | wrong: the record has no primary data behind it |

**Recommendation** (for the coordinator, not done here): replace the eight numbers with Table 1 values and their real
conditions, or withdraw the three records; drop "trapping", "Mark-Houwink" and "melanoidin" labels; key the paper as
`zha2019` in `data/keys/papers.yml` with this dossier.

**Outcome (2026-11-02, owner: "do whatever is best long term").** The three records were WITHDRAWN (deleted from
`process_state_calibrations.json`, `computational_priors.json` and `benchmark_intake_registry.json`), not corrected:
their premise (conjugation raises trapping pressure) is the opposite of what the paper proposes, and the real
Table 1 values live in the ledger entry `zha2019_tab1_sec_malls_pph_gum_arabic_conjugate_molar_mass`. "Zha et al.
(2019)" and "J. Agric. Food Chem. 2019 (Ref. 21)" were added to `closed_citations` so the research queue does not
re-encode the corpus values. Precedent: the Hidalgo & Zamora 2004 prior, deleted 2026-10-27 after its dossier.

## 7. Proposed ledger entries (YAML, not written into the ledgers)

Ledger choice: free amines, A304/A420 and colour in **trunk.yml** (time-resolved browning and amine loss; the
aroma_other header also lists protein "lysine loss", hence the also_relevant flag); SEC-MALLS and the non-lipid volatile
rows in **aroma_other.yml**; lipid volatiles and emulsion oxidation in **lipid.yml** (also_relevant: matrix for the
retention question). None is simulable: the engine has no protein-polysaccharide dry-state system, no a_w state and no
open-vessel loss term.

```yaml
# trunk.yml
  - id: zha2019_fig3b_free_amino_groups_pph_gum_arabic_dry_state
    dossier: zha2019_extraction.md
    paper: "Zha, Yang, Rao & Chen 2019 J Agric Food Chem 67:10195"
    lab: "North Dakota State Univ. (Chen B / Rao J)"
    system: "pea protein hydrolysate (PURIS 870H, DH 12 %) alone, and PPH + gum arabic 1:4 w/w; hydrated, pH 7.0, freeze-dried; 5 g powder"
    medium: "dry state, open 60 mL jar in desiccator over saturated KBr (79 % RH set point), air"
    temperature_c: [60]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [0, 1440, 4320, 7200]
    species: [free amino groups]
    kind: time_series
    location: "Fig. 3B (dossier sec. 4c digitisation, mmol/g protein); text -26.5 % / -9.8 % at 3 d"
    status: digitised_in_dossier
    quant_class: unknown
    simulable: false
    condition_quote: "maintained at 79% by saturated KBr solution and 60 °C"
    flags: [assay_not_named_in_paper, unit_implausible_high, amine_loss_exceeds_ga_reducing_ends, aw_set_point_not_measured, open_vessel, low_moisture_no_aw_term, "also_relevant: aroma_other"]
    note: "Directional only: GA ~3x the amine loss of heated PPH, mostly in day 1. ~6.2 mmol/g protein at t0 is ~4x the DH-12 % ceiling; GA-attributable loss ~180x the GA reducing ends (dossier sec. 1)."
  - id: zha2019_fig3a_a304_a420_pph_gum_arabic_dry_state
    dossier: zha2019_extraction.md
    paper: "Zha, Yang, Rao & Chen 2019 J Agric Food Chem 67:10195"
    lab: "North Dakota State Univ. (Chen B / Rao J)"
    system: "PPH alone and PPH + gum arabic 1:4 w/w, freeze-dried"
    medium: "dry state, open jar, 79 % RH (saturated KBr), air; absorbance of a redissolved solution of unstated concentration"
    temperature_c: [60]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [0, 1440, 4320, 7200]
    species: [A304, A420]
    kind: time_series
    location: "Fig. 3A (dossier sec. 4d digitisation, two y-axes)"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "measuring the absorbance of glyco-PPH solution at 304 and 420 nm"
    flags: [sample_concentration_not_stated, a304_not_a294, aw_set_point_not_measured, open_vessel, low_moisture_no_aw_term]
    note: "A304 (glyco) peaks at 3 d then falls; A420 rises monotonically to ~0.18 at 5 d vs ~0.045 for PPH alone."
  - id: zha2019_fig3c_colour_lab_pph_gum_arabic_dry_state
    dossier: zha2019_extraction.md
    paper: "Zha, Yang, Rao & Chen 2019 J Agric Food Chem 67:10195"
    lab: "North Dakota State Univ. (Chen B / Rao J)"
    system: "PPH alone and PPH + gum arabic 1:4 w/w, freeze-dried powders"
    medium: "dry state, open jar, 79 % RH, air; Minolta CR-310 on powder"
    temperature_c: [60]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [0, 1440, 4320, 7200]
    species: [L*, b*, a*]
    kind: time_series
    location: "Fig. 3C bubble chart (L* vs b*; a* as bubble size, not readable); dossier sec. 4e"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: ""
    flags: [a_star_not_readable, ga_dilutes_colour_at_t0, aw_set_point_not_measured]
    note: "L* 83.6 -> 72.4 (PPH+GA) vs 82.8 -> 78.1 (PPH alone) over 5 d."

# aroma_other.yml
  - id: zha2019_tab1_sec_malls_pph_gum_arabic_conjugate_molar_mass
    dossier: zha2019_extraction.md
    paper: "Zha, Yang, Rao & Chen 2019 J Agric Food Chem 67:10195"
    lab: "North Dakota State Univ. (Chen B / Rao J)"
    system: "PPI, PPH, gum arabic, PPH+GA 1:4 mixture (0 d) and glyco-PPH after 1 / 3 / 5 d at 60 C / 79 % RH"
    medium: "10 mM PBS pH 7.0, soluble fraction after 2000 g and 0.45 um filtration; SEC-UV-RI-MALLS, 22 C"
    temperature_c: [60]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [0, 1440, 4320, 7200]
    species: [Mw, Mn, Rg, Mw/Mn, mass fraction, N gum arabic per PPH]
    kind: endpoint
    location: "Table 1 (16 fraction rows; dossier sec. 4a); N ~1.2 in text and abstract"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "The gum arabic refractive index increment (dn/dc) was set at 0.141 mL/g"
    flags: [single_dndc_for_protein_and_polysaccharide, insoluble_fraction_discarded, n_not_a_graft_ratio, n_exceeds_available_ga, mass_recovery_column_over_100, no_viscometry_no_mark_houwink]
    note: "The only real numbers behind the three repo records jafc_2019_ref21_pea_gum_arabic_architecture_* (whose values are not in the paper; dossier sec. 6). Fraction-1 Mw 2.48 MDa (0 d) -> 4.17 (1 d) -> 3.66 -> 3.36 MDa."
  - id: zha2019_tab2_nonlipid_beany_markers_dry_incubation
    dossier: zha2019_extraction.md
    paper: "Zha, Yang, Rao & Chen 2019 J Agric Food Chem 67:10195"
    lab: "North Dakota State Univ. (Chen B / Rao J)"
    system: "PPH, gum arabic, PPH+GA 1:4 (0 d) and glyco-PPH 1 / 3 / 5 d (60 C, 79 % RH, open jar)"
    medium: "powder, HS-SPME-GC-MS (parameters in Zha 2019 Food Hydrocoll ref. 16), external standard calibration"
    temperature_c: [60]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [0, 1440, 4320, 7200]
    species: [3-methyl-1-butanol, acetophenone]
    kind: time_series
    location: "Table 2 rows 2 and 6 (dossier sec. 4b)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: ""
    flags: [unit_footnote_contradicts_text, basis_per_sample_or_per_pph_unclear, no_heated_pph_alone_control, open_vessel, "also_relevant: matrix", "complements: zha2019_tab2_lipid_beany_markers_dry_incubation"]
    note: "Alcohol and aryl ketone that cannot react with amines fall by the same ~60-67 % in day 1 as hexanal: evidence for physical loss over Maillard consumption."
  - id: zha2019_si_tabs1_95_volatiles_glyco_pph
    dossier: zha2019_extraction.md
    paper: "Zha, Yang, Rao & Chen 2019 J Agric Food Chem 67:10195"
    lab: "North Dakota State Univ. (Chen B / Rao J)"
    system: "PPH, GA, PPH+GA (0 d), glyco-PPH 1 / 3 / 5 d"
    medium: "powder, HS-SPME-GC-MS, relative content normalised to E-2-octenal"
    temperature_c: [60]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [0, 1440, 4320, 7200]
    species: [alkylpyrazines, 2-methylthiophene, 2-ethylthiophene, benzaldehyde, furfural, 5-methylfurfural, "4,5-dihydro-5-alkyl-2(3H)-furanones"]
    kind: time_series
    location: "Supporting Information Table S1 (ACS, DOI 10.1021/acs.jafc.9b04099); main text names compounds only"
    status: supplement_missing
    quant_class: semi_quantitative
    simulable: false
    condition_quote: ""
    flags: [supplement_missing, internal_standard_normalised_only, open_vessel]
    note: "95 volatiles (7 pyrazines, 2 thiophenes); pyrazines/Strecker aldehydes rise with time per text. Needed only if a dry-state protein-polysaccharide pyrazine check is wanted."

# lipid.yml
  - id: zha2019_tab2_lipid_beany_markers_dry_incubation
    dossier: zha2019_extraction.md
    paper: "Zha, Yang, Rao & Chen 2019 J Agric Food Chem 67:10195"
    lab: "North Dakota State Univ. (Chen B / Rao J)"
    system: "PPH (~8 % lipid), gum arabic, PPH+GA 1:4 (0 d) and glyco-PPH 1 / 3 / 5 d (60 C, 79 % RH, open jar)"
    medium: "powder, HS-SPME-GC-MS, external standard calibration (LOD 0.008-0.025 mg/L, recovery 92-98 %)"
    temperature_c: [60]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [0, 1440, 4320, 7200]
    species: [hexanal, 2-pentylfuran, 1-octen-3-one, 1-octen-3-ol]
    kind: time_series
    location: "Table 2 rows 1, 3, 4, 5 (dossier sec. 4b)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "Units of milligrams per gram of dry weight"
    flags: [unit_footnote_contradicts_text, basis_per_sample_or_per_pph_unclear, no_heated_pph_alone_control, open_vessel, matrix_food, "also_relevant: matrix", "complements: zha2019_tab2_nonlipid_beany_markers_dry_incubation"]
    note: "Hexanal 3.17 -> 0.42 (ppm reading) over 5 d; loss mechanism (stripping vs retention vs consumption) not separable (dossier sec. 1 ans. 2)."
  - id: zha2019_fig6_emulsion_looh_hexanal_37c
    dossier: zha2019_extraction.md
    paper: "Zha, Yang, Rao & Chen 2019 J Agric Food Chem 67:10195"
    lab: "North Dakota State Univ. (Chen B / Rao J)"
    system: "2 wt % corn oil (PUFA ~57 %) emulsions stabilised by 0.20 wt % PPH, PPH+GA mixture, or glyco-PPH-1; 0.04 % NaN3"
    medium: "10 mM PBS pH 7.0, homogenised (5000/500 psi x3)"
    temperature_c: [37]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [0, 14400]
    species: [lipid hydroperoxides, hexanal]
    kind: time_series
    location: "Fig. 6A (LOOH, uM) and 6B (hexanal, ng/mL; LOD 7.89) (dossier sec. 4f digitisation)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "measuring the formation of LOOH and hexanal in emulsion during storage at 37 °C"
    flags: [x_axis_mislabelled_conjugation_time, emulsion_not_cooking, droplet_size_differs_between_arms]
    note: "Hexanal lag ~4 d (PPH), ~7 d (mixture), ~9 d (glyco-PPH-1). Droplet d43 20.7 / 5.08 / 0.75 um differ, so interface area is confounded with emulsifier."
```

## 8. Supplementary material

The paper has Supporting Information (header symbol and p. 10 "Associated Content"): **Table S1, SPME-GC-MS profiled
volatile compounds in PPH-GA conjugates at different conjugation times** (PDF, free at pubs.acs.org under DOI
10.1021/acs.jafc.9b04099). **Not on disk** (no `zha2019_supplementary*` in `data/articles/`). It is the only place
with numbers for the Maillard volatiles (7 pyrazines, 2 thiophenes, Strecker aldehydes, furanones) and the other 89
volatiles, as relative content against E-2-octenal. **Needed?** Only if a dry-state pyrazine / Strecker check is wanted
(P(it changes a model parameter) ~3 %; P(useful as a directional pyrazine-vs-time check) ~30 %). Nothing in the SI
bears on the repo's SEC / architecture records (its caption is volatiles only). Every method detail (HS-SPME
parameters, amine assay, A304/A420 sample preparation) is in refs 16 and 18 (Zha et al. 2019 Food Hydrocoll. 92:30;
Zha et al. 2019 Food Chem. 285:130), not on disk.
