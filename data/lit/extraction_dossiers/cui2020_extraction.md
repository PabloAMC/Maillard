# Cui, Kimmel, Zhou, Rao & Chen 2020 — EXTRACTION (file stem `cui2020`)
DOI: **10.1016/j.carbpol.2020.116546** (printed on p. 1 and in the PDF metadata). Carbohydrate Polymers **245** (2020)
116546 (article number; 11 pp.). Received 15 Mar 2020, accepted 29 May 2020, online 5 Jun 2020.
Crossref (2026-11-02, `api.crossref.org/works/10.1016/j.carbpol.2020.116546`): title "Combining solid dispersion-based
spray drying with cyclodextrin to improve the functionality and mitigate the beany odor of pea protein isolate";
authors Leqi Cui, Jennifer Kimmel, Leon Zhou, Jiajia Rao, Bingcan Chen; issued 2020-10; Carbohydrate Polymers 245,
116546. **Matches the PDF and the brief** (publisher page S0144861720307207).
[!] **Naming caveat: the same first author and lab published at least three "Cui et al. 2020" papers on pea isolate**:
this one (Carbohydr Polym 245:116546); Cui, Kimmel, Zhou, Rao & Chen 2020 J Agric Food Res 100032 (extraction pH and
cultivar vs aroma, untargeted metabolomics; cited here as the MS-processing and seed-source method); and Cui,
Bandillo, Wang, Ohm, Chen & Rao 2020 Food Hydrocoll 108:106008 (cultivar and extraction pH vs functionality). The
stem `cui2020` is this paper only.

### Lab yellow-pea isolate (alkaline extraction pH 9.0 / isoelectric pH 4.5, not defatted, 85 % protein, 2.0 % lipid) spray-dried at 130 / 70 °C with or without α-, β- or γ-cyclodextrin at CD:DM 5:95 or 10:90, against the same isolate physically mixed with CD after drying. Volatiles by HS-SPME-GC-MS (60 °C, 1 g powder + 3 g 20 % NaCl), untargeted, shown **only as a log2 colour heat map** (23 compounds, no numbers). A spike test (hexanal, 1-octen-3-ol, 1-octen-3-one at 0.167-8.35 ppm on 1 g starch or pure α/β/γ-CD) gives **raw SPME peak areas** (Fig. 5). **No absolute volatile levels, no internal standard, no binding constant, no heating series, no LOX, no sensory panel.**

**Read 2026-11-02 from `data/articles/cui2020.pdf`** (11 pp., Elsevier, born-digital). Prose from `pdftotext -layout`
(scratchpad `ingest/cui2020_work/cui2020.txt`). All 11 pages rendered at 200 dpi (`p-01..p-11.png`); **Table 1 (p. 3)
and Table 3 (p. 7) read by eye from the 200 dpi renders**; **Fig. 4B heat map (p. 8) read from a 300 dpi crop**
(`fig4B.png`); **Fig. 5A-C (p. 9) digitised from 400 dpi crops** (`fig5A.png`, `fig5B.png`, `fig5C.png`), read error
±0.03 × 10⁷ (panel B, axis 0-2.5 × 10⁷) and ±0.1 × 10⁷ (panels A, C, axis 0-8 × 10⁷). Tables 2 and 4 and Figs. 1-3, 6
(structure, functionality) were looked at, not transcribed. Chen 2022's reference list read by `pdftotext` of
`data/articles/Chen2022.pdf` (sec. 1, answer 5).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure
(read error stated) · **[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]**
derived by me, never printed · **[!]** flag.

**Laboratory.** Department of Plant Sciences, **North Dakota State University**, Fargo (Cui, Rao, Chen) with **Roquette
America Inc.**, Geneva IL (Kimmel, Zhou). **Corresponding author Bingcan Chen** (bingcan.chen@ndsu.edu); Jiajia Rao
co-senior (conceptualisation, supervision, funding). Funded by Roquette [M, acknowledgement]; two Roquette co-authors
[!] (supplier interest in a commercial pea-protein process). **Not in `data/keys/laboratories.yml`: NEW.** The group
already appears in the ledgers under **three free-text spellings** — "North Dakota State Univ. / USDA-ARS Fargo (Rao)"
(matrix.yml, `gao2020_tab2_...`), "North Dakota State University / USDA-ARS Fargo (Rao)" (lipid.yml, `gao2020_fig5_...`,
`gao2020_fig6_...`) and "North Dakota State (Rao J.)" (aroma_other.yml, `gao2020_fig5_..._methoxypyrazine...`) — with no
group entry, so any fit table using them would be refused today. Proposed id **`ndsu_rao_chen`** (entered 2026-11-02 as `ndsu_chen_rao`, string "North Dakota State Univ. (Chen B / Rao J)"), strings: the three
above plus "North Dakota State University (Chen B / Rao J)"; evidence: Gao 2020 (Rao corresponding; Cui and Chen
co-authors), this paper (Chen corresponding; Rao co-senior), shared HS-SPME method (Xu, Jin, Lan, Rao & Chen 2019) and
shared isolate protocol. Rao and Chen co-run one group (confidence ~90 %); merging is the conservative leave-one-out
choice. **Same-lab dossiers:** `gao2020_extraction.md` (same AE-IEP yellow-pea isolate, LOX and beany markers;
Cui is a co-author). `gao2020b`, `cui2022`, `cui2026` dossiers are other labs (checked headers: not NDSU).

## 1. What this paper answers for the model

| model need | what Cui 2020 gives | tag |
|---|---|---|
| baseline lipid-oxidation volatiles in a pea isolate (absolute, µg/kg) | **none**. 23 "representative" volatiles as a log2 colour heat map (Fig. 4B); no table of abundances | [F] colour only |
| which beany markers are present | hexanal, nonanal, 2-octenal, decanal, 1-octen-3-ol, 1-octen-3-one, 2-pentylfuran, hexanoic acid, benzaldehyde, pentadecanal (+ 13 others) in Table 3 | [M] |
| lipid substrate in the isolate | **lipid 2.03 ± 0.04 % (Soxhlet), protein 84.97 ± 2.99 %, moisture 5.08 %** (Table 1, as-is powder) | [M] |
| CD effect on headspace of beany volatiles in PPI | heat map only; co-spray-dried PPI-CD lower than control for 21 of 23 compounds, physically mixed ≈ control | [F] / [Q] |
| inclusion / binding constants (K, CD-volatile) | **none** | — |
| partition / retention vs a reference matrix | Fig. 5: SPME peak area of 3 odorants spiked on pure CD vs on corn starch, 4 levels | [F] |
| PPI-matrix spike recovery (protein + CD) | Fig. S3, **SI only, not on disk** | — |
| heating / Maillard / T series | spray drying at one condition (inlet 130 °C, outlet 70 °C); **no T series**; ethylpyrazine higher in CD samples, Maillard **speculated** | [M] / [P] |
| LOX activity | **not measured** | — |
| sensory | **none**; "beany odour mitigation" is inferred from the instrumental heat map | [!] |

**Answers.**

1. **Volatiles measured, units, native PPI levels.** Untargeted HS-SPME-GC-MS; after MPP filtering 23 compounds were
   kept (Table 3): 3 alcohols, 7 aldehydes, 2 ketones, 3 acids, 2 alkenes, 1 ester, 2 heterocyclics, 3 unknowns.
   Among the brief's list, **hexanal, nonanal, 1-octen-3-ol and 2-pentylfuran are in; (E,E)-2,4-decadienal,
   (E)-2-nonenal, (E,E)-2,4-nonadienal, octanal, pentanal, 1-hexanol and 1-pentanol were run as authentic standards
   but did not survive the filter** [M, sec. 2.6.3] (i.e. not detected in all three replicates of any group, CV ≥ 25 %,
   or no ≥ 2-fold ANOVA-significant difference between groups — the paper does not say which). Standards served
   identification only; **quantification basis: none** — abundances are log2-normalised MS feature areas, shown as
   colours on a ±14.3 scale without a printed baseline. **No native level in any unit is printed** (confidence ~98 %).
   The only numeric native datum relevant to the lipid lane is the **isolate lipid content, 2.03 ± 0.04 % w/w as-is**
   (Gao 2020, same lab, same protocol class, printed 1.47 %). **Not usable as a lipid-lane starting-state check for
   volatiles** (confidence ~95 %); usable at most as a lipid-substrate figure for a lab AE-IEP pea isolate.
2. **CD reduction per volatile, ratio, and mechanism.** α-, β- and γ-CD at **CD:DM 5:95 ("-L") and 10:90 ("-H")**,
   i.e. CD = 5 or 10 % of feed dry matter; ≈ 0.124 g CD per g protein at 10:90, ≈ 0.10-0.13 mmol CD per g protein [Z,
   from Table 1 protein 84.97 % at 5.08 % moisture]. **No per-volatile reduction factor is printed.** The heat map
   (Fig. 4B, sec. 4c) puts the four controls at the top of the colour scale (+14.3, saturated red) and every co-spray-
   dried group near 0 for hexanal, nonanal, 2-pentylfuran, benzaldehyde, pentadecanal, valencene and α-farnesene; for
   2-octenal, decanal, acetoin, the acids, 3-methyl-1-butanol and 1-octen-3-ol the controls sit near 0 and some CD
   groups at −8 to −14 (blue). A ±14.3 log2 span is a 2 × 10⁴-fold ratio, not credible for 5-10 % CD when pure CD
   (Fig. 5) lowers the same odorants only 2-100 ×; the saturated cells almost certainly encode **non-detection with MPP
   missing-value imputation** (confidence ~70 %). So the heat map reads as "hexanal etc. present in controls, absent or
   near-absent in co-spray-dried PPI-CD", nothing finer.
   **Trapping (partition) vs process loss: the paper cannot distinguish them, and its own control argues against the
   authors' partition reading.** Every vial was the same slurry (1 g powder + 3 g 20 % NaCl, 60 °C, 60 min), so the
   physically mixed samples contained the same CD at the same 10:90 ratio in the vial as the co-spray-dried ones, yet
   grouped with the control in PCA and HCA (Fig. 4A/B) [Q/F]. A reversible inclusion equilibrium in solution would act
   identically in both. The extra drop in co-spray-dried samples must therefore come from a **lower total volatile
   content leaving the dryer** (stripped to the exhaust, or held in a form that does not re-equilibrate in a 60-min
   60 °C slurry) — i.e. a **process effect, not a matrix-partition effect** (confidence ~75 %). The authors call it CD
   "entrapping aroma compounds during spray drying" [V, abstract] but report no mass balance (no solvent extraction,
   no CD dissolution/release test of total content). Also against strong in-vial trapping: at 0.167 ppm spikes, the
   physical mix lowered headspace vs PPI-control for 1-octen-3-ol and 1-octen-3-one **but not for hexanal** [Q, text
   on Fig. S3]. At ~30 mM CD in the vial slurry [Z], any effective CD-hexanal K above ~30 M⁻¹ would halve free hexanal;
   the null therefore bounds the effective K in this protein slurry low (confidence ~55 %, since the heat map's
   saturation could hide a 2-fold change and protein competition is uncharacterised). **No binding or inclusion
   constant anywhere** (confidence ~99 %).
3. **Heating, Maillard, LOX, T series: none usable.** One spray-drying condition (inlet 130 °C, outlet 70 °C; residence
   time not stated). **Ethylpyrazine was higher in all CD co-spray-dried groups than in controls** (heat map; the
   authors' only "increase") and they "speculated" a CD-facilitated Maillard reaction during drying [V/P]; untested,
   no sugar/amine data, and pyrazines are also native to peas [C, Roland 2017]. No LOX assay, no inactivation step; the
   60 °C slurry during SPME could itself form LOX volatiles (same lab measured residual LOX in a finished isolate, Gao
   2020) [!]. No temperature or time series (confidence ~99 %).
4. **Ledger usability.** (a) **matrix.yml**: Fig. 5 is a semi-quantitative, three-odorant × four-level comparison of
   headspace peak area over pure α/β/γ-CD vs corn starch (dry or slurried powder — ambiguous, sec. 5). It is a
   relative suppression ratio against another solid, not a partition coefficient (no empty-vial or water reference,
   no headspace calibration). Proposed as `digitised_in_dossier`, not simulable. Confidence it can set any constant in
   `parameters_matrix.py`: ~3 %. (b) **lipid.yml**: the Fig. 4B heat map is figure-only colour data on volatiles in a
   pea isolate; indexed for completeness (as `gao2020_fig5_beany_markers_peak_areas` was), not usable as a baseline.
   Isolate lipid content 2.03 % goes in its note. P(either entry changes a shipped number) ~2 %.
5. **Is this the "Cui 2020 method" of `chen2022_extraction.md`? No** (confidence ~97 %). Chen 2022's reference list
   gives **Cui, Bandillo, Wang, Ohm, Chen & Rao 2020, Food Hydrocolloids 108:106008** (pdftotext of `Chen2022.pdf`).
   The protocols also differ: Chen describes hexane defatting, pH 8.5, 2 h, 25 °C, freeze-drying; this paper uses
   no defatting, pH 9.0, 1 h, 23 °C, spray-drying. Same NDSU lab, different paper.

## 2. Conditions

- **Raw material** [M]: yellow pea seeds, NDSU pulse crop breeding programme, crop year 2018; **cultivar not named**
  although chosen "to avoid differences caused by different cultivars" [V] [!]. Dehulled (roll stands), milled (Retsch
  ZM 200) with vacuum cooling to avoid heating the flour. **Not defatted.**
- **Isolate (PPI)** [M]: flour:water 1:15 (w/v), pH 9.0 with 1.0 M NaOH, stirred 1 h at 23 °C; 5524 × g 20 min 4 °C;
  supernatant filtered (Whatman 1), **pH 4.5** (acid not named), 5524 × g 10 min; pellets of several batches pooled;
  pellet DM 19.61 ± 0.24 %; resuspended, **pH 7.0**, feed DM adjusted to **13.0 %** (Eq. 1). No curd washing stated.
- **CD** [M]: Wacker CAVAMAX W6 (α), W7 (β), W8 (γ), ≥ 98 %. Added to the resuspended feed at CD:DM 5:95 or 10:90 (Eqs.
  2-4), DM kept at 13.0 %.
- **Spray drying** [M]: Büchi B-290 mini; pump 0.3 L/h, aspirator "40 L/h" [!: B-290 aspirator is normally given in
  m³/h; unit likely misprinted], **inlet 130 °C, outlet 70 °C**. Residence time, feed temperature, atomising gas not
  stated. Single drying condition; batches "repeated at least twice".
- **Physical mixtures** [M]: spray-dried PPI-control powder + α/β/γ-CD weighed and mixed at 10:90 (aroma, FTIR,
  functionality) or 5:95 (functionality only).
- **Groups for aroma** (10 × triplicate = 30 files) [M]: control; αCD-H/-L; βCD-H/-L; γCD-H/-L; αCD-, βCD-, γCD-Mix (10:90).
- **Storage** of powders between drying and analysis: **not stated** [!].
- **Proximate** [M]: Kjeldahl protein (factor not stated), oven moisture, Soxhlet lipid, combustion ash, carbohydrate
  by difference (Table 1, sec. 4a).
- **Spike test (Fig. 5 / Fig. S3)** [M]: matrices corn starch, pure α-, β-, γ-CD (set 1); PPI-control, physical mix,
  co-spray-dried PPI-CD (set 2). Matrix frozen 2 h; 1 g pre-cooled matrix in a 20 mL vial; cold standards of
  1-octen-3-ol, hexanal and 1-octen-3-one added to a final 0.167, 0.835, 1.67, 8.35 ppm (basis — per g matrix, per
  vial contents — and solvent **not stated**) [!]; vortexed, sonicated 30 min, "equilibration for three days prior to
  analysis" [V] (**temperature not stated**); then HS-SPME "using the same procedure as described above" [V] — whether
  the 3 g 20 % NaCl was added is **ambiguous** (~60 % that it was) [!].
- **pH** in the SPME vials: not measured. **O2**: air headspace (20 mL vial). **Water activity**: not measured.

## 3. Analytical methods

- **HS-SPME** [M]: 1 g powder + 3 g 20 % NaCl in 20 mL vial (PTFE/silicone septum); 10 min incubation at 60 °C, 250 rpm
  (CTC PAL RSI 120); **DVB/CAR/PDMS** fibre 50 min at 60 °C; desorption 3 min at 250 °C.
- **GC-MS** [M]: Agilent 7890B / 5977A, ZB-Wax 60 m × 0.25 mm × 0.25 µm; 40 °C → 45 °C/min → 85 °C → 9 °C/min → 200 °C
  → 45 °C/min → 250 °C (12 min); He 2 mL/min; EI 70 eV, source 230 °C. Method per Xu, Jin, Lan, Rao & Chen 2019 [C].
- **Data processing** [M]: MassHunter deconvolution; Mass Profiler Professional B.14.5 alignment (RT window 0.15 min,
  spectral match factor 0.3); keep features present in 100 % of replicates of at least one group, CV < 25 % on log2
  abundances within group, one-way ANOVA p < 0.05 and fold change ≥ 2 → 23 features. PCA, HCA, heat map in MPP.
  Normalisation/baselining of the heat map (to median? to control?) **not stated**; missing-value handling **not
  stated** [!].
- **Identification** [M]: NIST 17 match > 80, RI vs NIST, and authentic standards for 7 of 23 (superscript a in Table
  3: 3-methyl-1-butanol, 1-octen-3-ol, hexanal, nonanal, 2-octenal, 1-octen-3-one, 2-pentylfuran). Others are
  library/RI-tentative.
- **Quantification basis: none.** No internal standard, no response factors, no calibration curves in matrix, no LOD.
  Fig. 5 is raw peak area (axis label "Peak aera" [sic]); its error bars are undefined (SD presumed, n = 3?).
- **Statistics** [M]: "repeated at least twice", triplicate samples; ANOVA + Tukey HSD (SPSS 25); ± in Tables 1, 2, 4
  undefined (presumably SD).

## 4. Data

### 4a. Table 1 (p. 3): proximate composition of spray-dried powders, % (as printed; read by eye at 200 dpi) [M]

| sample | CD:DM | protein | carbohydrate | moisture | lipid | crude ash |
|---|---|---|---|---|---|---|
| PPI-control | 0:100 | 84.97 ± 2.99 | 3.24 ± 1.87 | 5.08 ± 0.03 | 2.03 ± 0.04 | 4.69 ± 0.05 |
| PPI-αCD | 10:90 | 78.64 ± 1.06 | 9.58 ± 1.02 | 5.04 ± 0.01 | 2.00 ± 0.04 | 4.74 ± 0.02 |
| PPI-αCD-L | 5:95 | 83.55 ± 1.09 | 4.71 ± 1.15 | 5.07 ± 0.04 | 2.02 ± 0.07 | 4.66 ± 0.04 |
| PPI-βCD | 10:90 | 78.71 ± 0.78 | 9.56 ± 0.70 | 5.07 ± 0.02 | 1.99 ± 0.09 | 4.67 ± 0.00 |
| PPI-βCD-L | 5:95 | 82.94 ± 0.88 | 5.34 ± 0.59 | 5.08 ± 0.05 | 1.98 ± 0.08 | 4.67 ± 0.14 |
| PPI-γCD | 10:90 | 79.15 ± 0.47 | 9.06 ± 0.47 | 5.06 ± 0.00 | 2.02 ± 0.07 | 4.72 ± 0.08 |
| PPI-γCD-L | 5:95 | 83.09 ± 0.46 | 5.20 ± 0.63 | 5.06 ± 0.01 | 1.96 ± 0.05 | 4.68 ± 0.12 |

Dry-basis control [Z]: protein 89.5 %, lipid 2.14 %. Carbohydrate rows close by difference (e.g. 100 − 78.64 − 5.04 −
2.00 − 4.74 = 9.58) [Z].

### 4b. Table 3 (p. 7): the 23 representative compounds [M] (RI calculated / reference; a = authentic standard)

Alcohols: 3-methyl-1-butanol^a 1226/1221; 1-octen-3-ol^a 1459/1459; pentadecanol 2296/2287.
Aldehydes: 3-methylbutanal 910/913; hexanal^a 1057/1068; nonanal^a 1409/1408; 2-octenal^a 1456/1451 (geometry not
stated); decanal 1517/1516; benzaldehyde 1566/1557; pentadecanal 2050/2054.
Ketones: acetoin 1311/1309; 1-octen-3-one^a 1345/1331. Acids: acetic 1467/1466; 3-methylbutanoic 1687/1688; hexanoic
1868/1869. Alkenes: valencene 1760/1764; α-farnesene 1766/1756. Ester: ethyl isohexanoate 1186/1193. Heterocyclics:
2-pentylfuran^a 1232/1231; ethylpyrazine 1379/1371. Unknowns RI 806, 1399, 1444 (heat-map labels 78.0@7.894714,
57.0@11.84125, 175.0@12.390916; the mapping of RI to label is not printed).

### 4c. Fig. 4B heat map (p. 8), read from a 300 dpi crop [F] (colour classes; scale −14.3 … 0 … +14.3, log2)

Column order (HCA): α-farnesene, 78.0@7.89, valencene, 1-octen-3-one, benzaldehyde, pentadecanal, 2-pentylfuran,
nonanal, hexanal | acetoin, hexanoic acid, acetic acid, 3-methylbutanoic acid, 2-octenal, decanal, 3-methyl-1-butanol,
57.0@11.84, pentadecanol, ethylpyrazine, 1-octen-3-ol, 175.0@12.39, ethyl isohexanoate, 3-methylbutanal.
Row clusters: {βCD-Mix, αCD-Mix, Control, γCD-Mix} vs {γCD-H, γCD-L} and {βCD-H, βCD-L, αCD-H, αCD-L}.

| compound | Control + 3 Mix | co-spray-dried CD groups |
|---|---|---|
| hexanal, nonanal, 2-pentylfuran, benzaldehyde, pentadecanal, valencene, α-farnesene, 78.0@7.89 | +14.3 (saturated red; αCD-Mix pentadecanal ≈ +10) | ≈ 0 (pale; hexanal/nonanal slightly negative, ≈ −1 to −3, in γ and β rows); benzaldehyde ≈ −14 in γCD-H and βCD-H; 2-pentylfuran ≈ −7 in βCD-L |
| 1-octen-3-one | Control +14.3; α/β-Mix ≈ +10 | ≈ −8 in all except **αCD-H ≈ +10** |
| acetoin, hexanoic acid, acetic acid | ≈ 0 | ≈ −14 in α and β groups (hexanoic −7 in βCD-H); ≈ 0 in γ groups |
| 2-octenal, decanal, 3-methyl-1-butanol, 57.0@11.84, pentadecanol | ≈ 0 (decanal: βCD-Mix +14, Control/γCD-Mix ≈ +8) | ≈ −14 in γ-H/γ-L and several α cells; β mostly ≈ 0 |
| 1-octen-3-ol | ≈ 0 | βCD-H ≈ −14; βCD-L, γCD-L ≈ −7; others ≈ 0 |
| ethylpyrazine | ≈ −2 (pale green) | ≈ 0 to +4 (pale orange): **the only compound higher with CD** |
| 3-methylbutanal, ethyl isohexanoate, 175.0@12.39 | ≈ 0 | ≈ −14 in αCD-L; −7 in several others; ethyl isohexanoate ≈ +4 in αCD-H, β rows |

Read error: one colour class (≈ ±3 log2 units); the scale saturates at both ends, so cells at ±14.3 have no upper
bound on magnitude.

### 4d. Fig. 5 (p. 9): headspace SPME peak area (× 10⁷) of spiked odorants on 1 g matrix [F]

Hexanal (panel B, read ±0.03):

| spike (ppm) | starch | αCD | βCD | γCD |
|---|---|---|---|---|
| 0.167 | 0.19 | 0.05 | 0.10 | 0.09 |
| 0.835 | 0.62 (bar ±0.21) | 0.05 | 0.26 | 0.28 |
| 1.67 | 1.01 (±0.17) | 0.38 | 0.55 | 0.55 |
| 8.35 | 2.25 (no bar) | 0.45 (±0.16) | 1.14 (±0.31) | 1.05 (±0.36) |

1-octen-3-ol (panel A, read ±0.1):

| spike (ppm) | starch | αCD | βCD | γCD |
|---|---|---|---|---|
| 0.167 | 0.40 | 0.0 | 0.0 | 0.1 |
| 0.835 | 1.6 (bar 1.25-2.07) | 0.0 | 0.2 | 0.6 |
| 1.67 | 2.9 (2.05-3.82) | 0.04 | 0.3 | 1.0 |
| 8.35 | 7.5 (6.4-8.55) | 0.7 | 1.6 | 1.75 |

1-octen-3-one (panel C, read ±0.1):

| spike (ppm) | starch | αCD | βCD | γCD |
|---|---|---|---|---|
| 0.167 | 0.73 | 0.39 | 0.17 | 0.33 |
| 0.835 | 1.63 (bar 0.95-2.3) | 0.08 | 0.39 | 0.37 |
| 1.67 | 2.23 (1.95-2.5) | 0.32 | 0.17 | 0.58 |
| 8.35 | 6.2 (4.85-7.5) | 0.58 | 1.62 | 2.10 |

Derived CD/starch peak-area ratios [Z] (≥ 0.835 ppm, where the starch signal is ≥ 10 × read error):
hexanal α 0.08 / 0.38 / 0.20, β 0.42 / 0.54 / 0.51, γ 0.45 / 0.54 / 0.47 (0.835 / 1.67 / 8.35 ppm);
1-octen-3-ol α ≤ 0.06 / 0.01 / 0.10, β 0.13 / 0.11 / 0.21, γ 0.38 / 0.35 / 0.23;
1-octen-3-one α 0.05 / 0.14 / 0.09, β 0.24 / 0.08 / 0.26, γ 0.23 / 0.26 / 0.34.
Ordering: **α-CD suppresses most for all three odorants; β ≈ γ for hexanal** (~2 ×), α ~5 × at 8.35 ppm.
The starch curves are **concave, not linear** as the text says [Q vs F]: hexanal slope 0-1.67 ppm ≈ 0.6 × 10⁷/ppm,
1.67-8.35 ppm ≈ 0.19 × 10⁷/ppm [Z] (fibre competition on DVB/CAR/PDMS is a candidate cause).

### 4e. Not transcribed (functional/structural)

Table 2 (FTIR secondary structure %), Table 4 (WBC, OBC, emulsion particle size, ESI, foaming), Fig. 1 (SEM), Fig. 2
(FTIR), Fig. 3 (ζ-potential pH 3-8), Fig. 6 (solubility pH 3-8). No model use.

## 5. Flags and internal consistency checks

1. [!] **No absolute or semi-quantitative volatile numbers**; the aroma result is a colour heat map with an
   unstated baseline and saturated scale. Saturation at ±14.3 log2 (2 × 10⁴ ×) is implausible as a real ratio for
   5-10 % CD; most likely missing-value imputation (~70 %). A 2-fold difference between control and physical mix
   would be invisible in the saturated block.
2. [!] **The physical-mix control undercuts the partition interpretation** (sec. 1, answer 2): same CD dose in the vial,
   no effect; the co-spray-dried effect is a content difference with no mass balance. Authors' "CD entrapping aroma
   compounds during spray drying" is [P].
3. [!] **Table 1 dilution does not add up** [Z]. Adding 10 % CD should dilute protein, lipid and ash by 10 %: expected
   protein 76.5, lipid 1.83, ash 4.22 at 10:90. Printed protein 78.6-79.2 (≈ 75 % of the expected drop), **lipid
   1.99-2.02 and ash 4.67-4.74, i.e. no dilution at all**; carbohydrate by difference rises 5.8-6.3 points, not 10.
   Either less CD reached the powder than nominal, or lipid/ash were not measured per sample. Undermines "CD:DM
   10:90" as the realised dose (confidence ~65 % that the realised CD content is below nominal or the table is
   partly not per-sample).
4. [!] **1-octen-3-one in αCD-H is high (≈ +10)**, contradicting "all CD treatment … lighter color" [V, p. 7]; the
   recommended αCD at 5:95 is the low-dose arm, and αCD-L is indeed lower — dose ordering inverted for this compound.
5. [!] Text says controls are "a lighter color (blue)" for ethylpyrazine; the heat map shows them pale green (≈ −2),
   not blue — direction right, wording loose.
6. [!] Fig. 5 spike basis ("ppm") and solvent not stated; 3-day equilibration temperature not stated; whether NaCl
   solution was added before SPME ambiguous. If 3 g water was added, 1 g β-CD (aqueous solubility ≈ 18.5 g/L at 25 °C
   vs α ≈ 145, γ ≈ 232 g/L, standard values) would be mostly undissolved, so the α/β/γ comparison mixes solid-surface
   and solution complexation [Z].
7. [!] No empty-vial / water reference in Fig. 5: starch is itself a sorbent, so CD/starch ratios are relative to an
   uncharacterised solid, not to free partition.
8. [!] Starch "increased linearly" [V] vs concave data (sec. 4d).
9. [!] LOX not inactivated; 60 min at 60 °C in a 25 % w/w slurry can form or destroy LOX volatiles during sampling
   (same caveat as `gao2020_fig5_beany_markers_peak_areas`).
10. [!] No sensory evaluation; "mitigate the beany odor" is instrumental only.
11. [!] Roquette funding and two Roquette co-authors.
12. [!] Cultivar unnamed; storage time unstated; aspirator unit likely misprinted (L/h vs m³/h).
13. [!] MPP CV filter applied to log2 abundances (a CV < 25 % on a log scale is a weak reproducibility filter).
14. Three "Cui et al. 2020" papers from the same lab (header) — citation hazard.

## 6. Existing records in the repo that cite this paper

Searched 2026-11-02 with `grep -rIn -i` over `data/`, `docs/`, `src/`, `scripts/`, `tasks/` for
`carbpol.2020.116546`, `S0144861720307207`, `cui2020`, `cui 2020`, `cui et al`, `cui, l`, `cyclodextrin`,
`carbpol`, `co-spray`, `beany odor of pea`: **no record cites this paper.** Hits were other papers: Cui et al. 2022
JAFC 70:11556 (`computational_priors.json`, `benchmark_intake_registry.json`, `papers.yml`), Xin/Cui 2026 (ledgers),
Guo 2020 (β-CD flavour carrier, `guo2020_extraction.md`, matrix.yml), Weenen 1992 (α-CD HPLC column), Deng 2021 /
Zhai 2023c (Jiangnan "Cui et al." ARP method).

- `data/lit/extraction_dossiers/chen2022_extraction.md` l. 40, "**PPI (§2.2, Cui 2020 method)**": **not this paper**
  (Chen 2022 cites Cui, Bandillo, Wang, Ohm, Chen & Rao 2020 Food Hydrocoll 108:106008). Not a numeric mismatch,
  but the bare "Cui 2020" is ambiguous now that `cui2020` names a different paper; suggest the dossier spell out
  the Food Hydrocoll reference.
- `data/research_corpus/raw/Gemini_deep_research_08.md` l. 125 (β-CD "as low as 1 % w/w … ~60 %" headspace reduction
  of hexanal and (E)-2-nonenal) cites a 2025 ResearchGate review (its ref. 4), not this paper; **nothing in Cui 2020
  supports those numbers** (no 1 % dose, no (E)-2-nonenal in the final list, no percentage printed).
- `data/research_corpus/raw/Gemini_deep_research_I.md` l. 214 ("cyclodextrin washing") — no citation to this paper.

## 7. Proposed ledger entries (NOT written into the ledgers)

Two entries. The heat map is indexed because the lipid ledger indexes the same lab's equivalent figure (Gao 2020);
it carries no usable number. Fig. 5 is the only digitisable data set. The Fig. S3 protein-matrix spike set is SI
only (sec. 8) and gets no entry until the SI is on disk.

```yaml
  # ---- cui2020_extraction.md ----
  - id: cui2020_fig5_cyclodextrin_vs_starch_headspace_suppression
    dossier: cui2020_extraction.md
    paper: "Cui, Kimmel, Zhou, Rao & Chen 2020 Carbohydr Polym 245:116546"
    lab: "North Dakota State Univ. (Chen B / Rao J)"
    system: "hexanal, 1-octen-3-ol, 1-octen-3-one spiked at 0.167 / 0.835 / 1.67 / 8.35 ppm onto 1 g corn starch or pure alpha-, beta-, gamma-cyclodextrin"
    medium: "1 g matrix in 20 mL vial, equilibrated 3 days, then HS-SPME DVB/CAR/PDMS 60 C (3 g 20 % NaCl added: ambiguous)"
    temperature_c: [60]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [60]
    species: [hexanal, 1-octen-3-ol, 1-octen-3-one]
    kind: ratio
    location: "Fig. 5A-C (raw peak areas, 4 matrices x 4 spike levels; digitised in dossier sec. 4d with CD/starch ratios)"
    status: digitised_in_dossier
    quant_class: semi_quantitative
    simulable: false
    condition_quote: "followed by equilibration for three days prior to analysis"
    flags: [peak_area_only, no_free_partition_reference, spike_basis_not_stated, equilibration_temperature_not_stated, nacl_addition_ambiguous, beta_cd_solubility_limited, starch_curve_concave_not_linear, supplier_coauthors, "also_relevant: lipid"]
    note: "Relative headspace suppression over pure CD vs corn starch, not a partition or inclusion constant. alpha-CD strongest for all three (CD/starch 0.01-0.2 for 1-octen-3-ol/-one, 0.08-0.38 hexanal); beta ~ gamma ~0.5 for hexanal. No protein. The PPI-matrix spike set (PPI, physical mix, co-spray-dried PPI-CD) is Fig. S3, SI not on disk; text says physical-mix CD at 10 % of DM did not lower spiked hexanal headspace."
```

```yaml
  # ---- cui2020_extraction.md ----
  - id: cui2020_fig4_heatmap_pea_isolate_cyclodextrin_volatiles
    dossier: cui2020_extraction.md
    paper: "Cui, Kimmel, Zhou, Rao & Chen 2020 Carbohydr Polym 245:116546"
    lab: "North Dakota State Univ. (Chen B / Rao J)"
    system: "lab AE-IEP yellow pea isolate (pH 9.0 / pI 4.5, not defatted, 85 % protein, 2.03 % lipid) spray-dried 130/70 C alone or with alpha/beta/gamma-CD at CD:DM 5:95 or 10:90, vs physical mixtures 10:90"
    medium: "1 g powder + 3 g 20 % NaCl in 20 mL vial, HS-SPME DVB/CAR/PDMS 60 C, 10 + 50 min"
    temperature_c: [60]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [60]
    species: [hexanal, nonanal, 2-octenal, decanal, 1-octen-3-ol, 1-octen-3-one, 2-pentylfuran, hexanoic acid, benzaldehyde, ethylpyrazine]
    kind: ranking
    location: "Fig. 4B heat map (log2 colour, +-14.3, 23 compounds x 10 groups; colour classes in dossier sec. 4c); Table 3 identities; Table 1 proximate"
    status: figure_only
    quant_class: other
    simulable: false
    condition_quote: "One-gram powder sample was added with 3 g of 20 % NaCl solution"
    flags: [peak_area_only, colour_scale_saturated_likely_imputation, heatmap_baseline_not_stated, physical_mix_control_contradicts_partition_reading, no_mass_balance, lox_not_inactivated, matrix_food, cd_dilution_not_seen_in_table1, no_sensory, supplier_coauthors]
    note: "Qualitative only: hexanal, nonanal, 2-pentylfuran present in control and physical mixes, near-absent in co-spray-dried PPI-CD; ethylpyrazine higher with CD (Maillard speculated, untested). ph.initial is the spray-dryer feed pH. Isolate lipid 2.03 +- 0.04 % as-is (Soxhlet), protein 84.97 % (Table 1) -- the only lipid-lane number. No native volatile level in any unit."
```

## 8. Supplementary material

The paper has **Appendix A. Supplementary data** at the publisher (doi:10.1016/j.carbpol.2020.116546, ScienceDirect
"Supplementary material" of article S0144861720307207). Referenced items: **Fig. S1** (SEM 5000 ×), **Fig. S2** (FTIR
of β-/γ-CD samples), **Fig. S3** (spike test on PPI-control, physical mixture and co-spray-dried PPI-CD, with insets
at 0 ppm showing the **native peak areas of hexanal, 1-octen-3-ol and 1-octen-3-one** in each matrix), **Fig. S4**
(solubility of physical mixtures), **Table S1** (functionality at 5:95 and of physical mixtures). **Not on disk.**
Model-relevant: **Fig. S3 only** — it is the one place with (relative) native levels of three beany markers in the
isolate and the only protein-matrix ± CD spike data. Still peak areas without IS, and the authors say the protein
matrix gave no linear response. P(Fig. S3 yields a usable matrix or lipid-lane number) ~10 %; worth fetching only if a
CD-in-formulation question becomes active.
