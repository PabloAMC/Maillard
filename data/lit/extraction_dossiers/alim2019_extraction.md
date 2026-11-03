# Alim, Song, Liu, Zou, Zhang, Zhang & Raza 2019 — EXTRACTION (file stem `alim2019`)
LWT - Food Science and Technology **112**:108239, DOI 10.1016/j.lwt.2019.06.006. Crossref 2026-11-03: "Research of
beef-meaty aroma compounds from yeast extract using carbon module labeling (CAMOLA) technique", authors Aygul Alim,
Huanlu Song, Ye Liu, Tingting Zou, Yu Zhang, Songpei Zhang, Ali Raza; LWT vol. 112, article 108239, issued 2019-09.
**Matches the PDF and the brief.** Received 20 Feb 2019, online 7 Jun 2019. No naming caveat. Not to be confused with
the companion `alim2018` (JSFA 2018, thermally treated yeast extracts; queued, not read here).

### [13C5]-D-xylose : D-xylose 1 : 1 (50 + 50 µmol) + one synthetic yeast-extract peptide (100 µmol) ± L-cysteine 50 ± thiamine 50 µmol, "pH 6.0", 135 °C, 90 min, foil-sealed vial inside a Parr high-pressure reactor; HS-SPME-GC-O-MS end point. CAMOLA isotopologue % (Table 3) for 15 odorants in three labelled groups (V peptide-xylose, VI + Cys, VII + thiamine). Unlabelled volatiles (Table 2) for a yeast-extract pot and four peptide pots, internal-standard ratio, **no unit printed**. **The yeast extract itself was never labelled. No group with label + Cys + thiamine is reported (Group VIII run, no data). No kinetics, one T, one t, no SD, no [M−H]+ or natural-abundance correction.**

**Read 2026-11-03 from `data/articles/alim2019.pdf`** (8 pp., Elsevier, born-digital). Prose from
`pdftotext -layout` (scratchpad `ingest/alim2019_work/alim2019.txt`). All pages rendered at 150 dpi (`p-1.png` to
`p-8.png`). **Table 1 (p. 2), Table 2 (p. 6) and Table 3 (p. 7) read by eye from 300 dpi renders** (`hi-6.png`,
`hi-7.png`, Table 3 crop `t3crop.png`); the text layer agrees with the renders cell for cell for Tables 2-3 (empty
cells in Table 3 were checked on the crop). **Fig. 1 (ethylpyrazine spectra) and Fig. 2 (furfural spectra), p. 3,
digitised from 400 dpi crops** (`f400-3.png`, `ep_axis.png`, `furf_axis.png`), read error ±0.03 of the y-axis unit;
Fig. 3 (5-ethenyl-4-methylthiazole spectra, p. 4) from a 400 dpi crop (`thz_zoom.png`). Figs. 4-6 are pathway
schemes [P]. Sums, label shares and [M−H]+ re-allocations by hand [Z].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Beijing Advanced Innovation Center for Food Nutrition and Human Health / Beijing Engineering and
Technology Center of Food Additives, **Beijing Technology and Business University**; corresponding author **Huanlu
Song** (songhl@th.btbu.edu.cn). **Not in `data/keys/laboratories.yml`** (only `btbu_xie`, Xie Jianchun, whose
evidence line already says other BTBU groups are distinct). The ledgers already carry the free string
**"BTBU (Song H)"** / "BTBU (Song H.)" for `wang2012`, `xu2013` (sulfur.yml) and `raza2020` (aroma_other.yml), and
`xu2013_extraction.md` identifies Song's "Laboratory of Molecular Sensory Science" as distinct from the Xie group.
**NEW lab id proposed: `btbu_song`**, strings ["BTBU (Song H)", "BTBU (Song H.)"], evidence: Wang, Yang & Song 2012
Food Chem 131:280; Xu, Liu, Song et al. 2013 FRI 54:683; Raza, Song et al. 2020 Food Funct 11:8583; this paper — Song
corresponding in all four, same e-mail, same Parr-reactor + 2-methyl-3-heptanone-IS + 1:1 [13C5]-xylose design (Wang
2012, Xu 2013). Same-lab dossiers: `wang2012_extraction.md`, `xu2013_extraction.md`, `raza2020_extraction.md`. For
leave-one-lab-out, this paper is **not independent** of Wang 2012 / Xu 2013 (P ~0.95 same group and method).

## 1. What this paper answers for the model

| model need | what Alim 2019 gives | tag |
|---|---|---|
| MFT carbon origin, pentose + Cys, **no thiamine** (E6) | Group VI: M 25 / M+4 34 / M+5 41 %; no M+1..M+3 | [M] |
| MFT carbon origin with thiamine present (E6: thiamine share) | **not given**: labelled Group VIII (all four reactants) is in Table 1 but has no Table 3 rows | [!] |
| FFT, 2-methylthiophene carbon origin | Group VI: FFT 52 / 48 (M / M+5); 2-methylthiophene 42 / 48 (sum 90) | [M] |
| furfural, 2-methylfuran, 2-furanmethanol, 2-acetylfuran | Group V / VII isotopologue rows | [M] |
| pyrazine carbon modules (E13/E14) | Group V: methylpyrazine, 2,5-/2,6-DMP, 2-ethyl-5-methylpyrazine, trimethylpyrazine; Group VI ethylpyrazine | [M] |
| pyrazine nitrogen origin | **none** (no 15N; peptides are the only N source besides Cys/thiamine) | — |
| thiazole from thiamine | Group VII: 5-ethenyl-4-methylthiazole 100 % M | [M] |
| [M−H]+ / natural-abundance correction | **neither**; Figs. 1-2 show the [M−H]+ ion was counted as an isotopologue | [!] |
| yeast-extract composition (free AA, Cys, sugars, nucleotides, thiamine) | **not printed**; YE recipe not printed either (maybe SI Table 1) | [!] |
| amounts of odorants | Table 2, IS ratio (RF = 1), unit not printed, no SD, internal inconsistencies | [M] / [!] |
| kinetics, T or t series, pH after heating | **none** | — |

**Answers to the brief's questions.**

1. **Exact system** (confidence ~90 % on what is printed; much is not printed).
   - **Sugar is D-xylose, not ribose or glucose**: [13C5]-xylose (99 %, Isotec / Sigma-Aldrich) **1 : 1 with
     unlabelled xylose**, 50 + 50 µmol (Table 1) [M].
   - **The labelled pots are synthetic-peptide model systems, not yeast extract.** "Peptide 100 µmol" of
     Pro-Glu-Thr, Glu-Asn or Arg-Leu (95 %, Shanghai Jier); Table 1 does not say which peptide, or whether 100 µmol
     is one peptide or the mixture [!]. Cys 50 µmol (VI, VIII), thiamine 50 µmol (VII, VIII) [M].
   - **Yeast-extract pot**: YE paste LA00 (Angel Yeast, Yichang) + xylose + thiamine + Cys, **amounts not printed**
     (perhaps in Supplementary Table 1, sec. 8); **no YE composition** (free amino acids, Cys, reducing sugars,
     nucleotides, thiamine: none printed). The YE pot was run **unlabelled only** (Table 2 column) [M].
   - 135 °C, 90 min, "pH 6.0" (how set and whether buffered not stated; chemicals list NaH2PO3 / Na2HPO3, i.e.
     phosphite as printed, almost certainly phosphate [!]); sample in a vial "sealed with aluminum foil" placed in a
     Parr high-pressure reactor; **water volume not printed** (8 mL taken for SPME, so ≥ 8 mL) → **no molar
     concentrations computable** [!]. Triplicate stated.
   - [!] **Labelled groups V-VIII carry 100 µmol total xylose; their unlabelled twins I-IV carry 50 µmol.** The
     labelled and unlabelled pots are therefore not the same chemistry (sugar : N doubled).
2. **Isotopologue distributions**: transcribed in full in sec. 4b (Table 3). **Corrections: none** (confidence
   ~90 %). Evidence: (a) no correction is described; the only formula is a "labelled rate" that doubles the labelled
   share for the 1:1 dilution; (b) **no M+1 cell is non-blank anywhere** although natural 13C alone gives ~5-8 % M+1
   for C5-C8 odorants [Z], so small peaks were dropped or ignored, method unstated; (c) **Fig. 2 (furfural) proves
   no [M−H]+ correction**: the reference spectrum shows m/z 95 ([M−H]+) ≈ m/z 96 (M+•); in the labelled sample the
   authors assign m/z 100 as "[M+ +4]" and draw a four-label structure, whereas 100 is the [M−H]+ of the M+5
   isotopologue (101) [F/Z]; (d) **Fig. 1 (ethylpyrazine, MW 108) labels m/z 107, the [M−H]+ base peak, as "[M+]"**
   and counts 111/112/113 as M+4/M+5/M+6 from 107 [F]. So the M+4 column of every furan row (2-methylfuran,
   furfural, 2-acetylfuran, MFT, 2-furanmethanol) is plausibly the [M−H]+ of M+5, not a C4-labelled molecule
   (P ~0.8 overall; ~0.95 for furfural, where Fig. 2 shows it; ~0.7 for MFT, whose spectrum is not shown). This is
   the Song group (BTBU), so it is evidence on field practice, **not** on whether Zhou 2024 (Jiangnan; EX-052)
   corrected.
3. **Carbon origin, calibrated.**
   - **MFT (Group VI, peptide + xylose + Cys, no thiamine): essentially all intact xylose C5** (P ~0.75 that
     ≥ 90 % of MFT carbon is an intact xylose C5 skeleton here). Only M, M+4, M+5 appear; no M+1..M+3, so no
     detectable recombination of labelled and unlabelled fragments (C2+C3, C1+C4) and no detectable non-sugar carbon
     source. Re-allocating M+4 as [M−H]+ of M+5 and adding the matching unseen [M−H]+ of M: unlabelled 38 % /
     labelled 62 % [Z]; furfural (which can only come from xylose) gives 47 / 53 by the same rule [Z]. A non-xylose
     carbon source would push MFT's unlabelled share *above* furfural's; it is below. **It cannot discriminate the
     engine's two pentose routes** (norfuraneol and 1,4-dideoxyosone both keep C5 intact) and **says nothing on the
     thiamine share** (P ~0.95 no information: no labelled thiamine + Cys group reported).
   - **FFT and 2-methylthiophene: intact xylose C5** (FFT 52 : 48, 2-MT 42 : 48 with 10 % missing) (P ~0.8).
   - **Pyrazines (Group V, peptide + xylose): ring and side-chain carbon from xylose fragments** (P ~0.8 that
     ≥ 80 % of pyrazine carbon is xylose-derived; peptide carbon would raise M above the binomial value and it does
     not). Best evidence: **methylpyrazine** M / M+2 / M+3 / M+5 = 18 / 23 / 25 / 34 vs 25 each for one C2 + one C3
     fragment drawn independently from a 1:1 pool; **2,5-DMP** 22 / 48 / 30 vs 25 / 50 / 25 for C3 + C3. The other
     rows fit no simple model (2,6-DMP 38 / 56 / 6; 2-ethyl-5-methylpyrazine sums to 111; trimethylpyrazine has M+4
     48 and no M+1, contrary to the authors' "2,5-DMP + C1"; ethylpyrazine is on the wrong mass axis, sec. 5).
     Nitrogen origin: not traced. **The text's ratios do not match its own table** (methylpyrazine "1:1:1:3",
     2,6-DMP "1:2:1", 2-ethyl-5-methylpyrazine and trimethylpyrazine "1:2:2:1") [!].
   - **Agreement with the repo's labelling corpus**: MFT and FFT here agree with Cerny 2007b's no-thiamine control
     (system F, MFT > 95 % xylose, `cerny2007b_t6_no_thiamine_control`), with Cerny & Davidek 2003 (MFT 93 % ribose
     in the norfuraneol competition, `cerny2003_t3_nf_ribose_route_shares`) and with Cerny & Briffod 2007 (FFT > 99 %
     xylose). FFT here (≈ 50 : 50, i.e. ~100 % xylose) **contradicts the same lab's earlier claims** (Xu 2013: FFT and
     MFT unlabelled; Wang 2012: FFT 80 : 20), both run with thiamine at ~1.2 mM; for MFT the earlier "unlabelled"
     results were thiamine pots and are not contradicted (no thiamine in Group VI). Pyrazine carbon from sugar C2/C3
     modules agrees qualitatively with the sugar-fragment picture of Weenen 1992 (`weenen1992_text_1_13c_glucose_pyrazine_label_incorporation`) and Zhou 2025b (`zhou2025b_table2_13c5_xylose_isotopomers`).
4. **Ledger entries**: three proposed (sec. 7); none fit-grade (no concentrations, unstated volume, peptide
   reactants the engine does not carry, uncorrected ion basis). Best use: a directional check that MFT / FFT /
   2-methylthiophene stay intact-C5 in a no-thiamine pentose + Cys pot, from a lab independent of Cerny.

## 2. Conditions

- **Reactants (Table 1, µmol)** [M]: I peptide 100 + xylose 50; II + Cys 50; III + thiamine 50; IV + Cys 50 + thiamine
  50; V-VIII = I-IV plus [13C5]-xylose 50. Peptides Pro-Glu-Thr, Glu-Asn, Arg-Leu (95 %, synthesised, "identified
  peptides" from YE per Alim 2018b / 2019 JSFA [C]); which peptide in which group: not stated [!].
- **YE pot** [M]: YE paste LA00 (Hubei Angel Yeast) with xylose, thiamine and Cys; amounts not printed.
- **Chemicals** [M]: thiamine, L-Cys, GSH, D-xylose, NaH2PO3, Na2HPO3 (Tianjin Guangfu); GSH is listed but used in no
  group [!].
- **pH** [M]: "pH 6.0", initial only, adjusting agent and buffer strength not stated; no pH after heating.
- **T / t** [M]: 135 °C, 90 min; heat-up not stated. **Vessel**: vial sealed with aluminium foil inside a Parr
  high-pressure reactor; headspace / O2 not stated. Heating medium not stated.
- **Water / volume**: not printed [!]. **Replicates**: "repeated three times"; no SD anywhere.

## 3. Analytical methods

- **Extraction** [M]: 8 mL sample + 1 µL IS (2-methyl-3-heptanone, 0.816 µg/µL) in a 20 mL vial, 55 °C 20 min, then
  CAR/DVB/PDMS 65 µm fibre 40 min; desorption 250 °C 5 min.
- **GC-O-MS** [M]: Agilent 7890A-7000B, DB-WAX 30 m × 0.25 mm × 0.25 µm, He; 40 °C 3 min → 230 °C (rate not printed)
  5 min → 260 °C 3 min; EI 70 eV, m/z 35-550; NIST 05 + LRI (± 20) + odour; three GC-O panellists.
- **Quantification** [M]: Ci = CIS × Ai / AIS (RF = 1). [!] The IS is given two ways: "500 ng/g" (sec. 2.5) vs 0.816
  µg into 8 mL ≈ 102 ng/mL (sec. 2.3) [Z]. Table 2 prints **no unit and no ±**.
- **CAMOLA** [M]: isotopologue % from EI spectra; ion basis (molecular ion vs [M−H]+, integrated ion chromatogram vs
  apex spectrum) not stated; "labelled rate" = Σ[M+n]L / ([M]UL + Σ[M+n]L) × 2 × 100 %, never reported per compound.
  Applying it to Table 3 gives > 100 % for most rows (MFT 150 %, furfural 142 %, 2-methylfuran 128 %) [Z] — the
  formula is incompatible with the table.

## 4. Data

### 4a. Table 2 (p. 6): "Analysis of beef-meaty aroma compounds", IS ratio, **unit not printed**, no SD [M]

| compound | RI | odour (as printed) | YE system | Group I | Group II | Group III | Group IV |
|---|---|---|---|---|---|---|---|
| 2-methylfuran | 837 | smoky | 45.11 | 7.83 | 15.16 | — | 3171.07 |
| 2-methylbutanal | 891 | cocoa | 97.43 | 2.22 | 1.11 | — | 292.12 |
| 2-methylthiophene | 1050 | roast beef | 12.77 | — | 4.91 | — | 1758.51 |
| hexanal | 1060 | fat | 3.21 | — | — | — | 1380.09 |
| 2-pentylfuran | 1149 | butter | 1.82 | — | — | — | 31.34 |
| methylpyrazine | 1186 | popcorn | 31.1 | 11.16 | 0.06 | 14.24 | 150.21 |
| octanal | 1248 | fat | 2.28 | — | — | — | 45.08 |
| 2-methyl-3-furanthiol | 1270 | meaty | 8.95 | — | 5.83 | — | 18.45 |
| 2,5-dimethylpyrazine | 1287 | roast beef | 24.95 | 35.49 | 0.82 | 29.69 | 123.51 |
| 2,6-dimethylpyrazine | 1306 | roast beef | 27.54 | 3.48 | 0.32 | 0.51 | 51.52 |
| ethylpyrazine | 1314 | peanut butter | 19.96 | — | 0.58 | 4.93 | 688.03 |
| 2-ethyl-5-methylpyrazine | 1309 | fruit, sweet | 35.28 | 27.63 | — | 23.85 | — |
| 2-ethyl-6-methylpyrazine | 1368 | nutty | 27.04 | — | — | — | — |
| nonanal | 1374 | fat | 20.69 | — | 0.05 | — | — |
| trimethylpyrazine | 1386 | roast | 6.38 | 8.53 | — | 0.26 | — |
| 2-furfurylthiol | 1408 | coffee | 41.12 | — | 7.71 | — | 249.12 |
| 2-ethyl-3,5-dimethylpyrazine | 1427 | roast | 26.36 | — | — | — | — |
| 3-ethyl-2,5-dimethylpyrazine | 1427 | potato | 29.57 | 1180.69 | 818.72 | 1514.04 | 33464.88 |
| furfural | 1436 | bread | 3.545 | — | — | — | — |
| 2,5-dimethylfuran | 1463 | beef, coffee | 52.46 | 20742697.3 | — | — | — |
| 2-acetylfuran | 1478 | balsamic | 8.47 | — | — | — | — |
| benzaldehyde | 1493 | burnt sugar | 73.27 | — | — | 8.96 | 340.23 |
| 5-ethenyl-4-methylthiazole | 1498 | nutty | 11.08 | — | — | 8.45 | — |
| 5-methyl-2-furaldehyde | 1535 | sweety | 12.29 | — | — | — | — |
| 2-furanmethanol | 1568 | cooked sugar | 67.22 | — | — | 62.88 | 15.73 |
| 2-methoxyphenol | 1824 | aromatic | 47.83 | — | — | — | — |
| 4-(2-furanyl)-3-buten-2-one | 1828 | nutty | 51.07 | 13.91 | — | 6.77 | — |
| 2-acetylpyrrole | 1875 | bread | 27.78 | — | — | — | — |
| indole | 2407 | burnt | 62.92 | — | — | — | — |

Groups: I peptide-xylose; II + Cys; III + thiamine; IV + Cys + thiamine (all unlabelled, 50 µmol xylose).

### 4b. Table 3 (p. 7): "Proportion of labeled carbon atoms in the molecule (%)" [M]

"—" = printed dash; "(blank)" = cell left empty on the page (checked on the 300 dpi crop). Row sums [Z].

| group | compound (C atoms) | M+ | +1 | +2 | +3 | +4 | +5 | +6 | +7 | +8 | sum [Z] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| V | 2-methylfuran (C5) | 36 | — | — | — | 19 | 45 | — | — | — | 100 |
| V | methylpyrazine (C5) | 18 | — | 23 | 25 | — | 34 | — | — | — | 100 |
| V | 2,5-dimethylpyrazine (C6) | 22 | — | — | 48 | — | — | 30 | — | — | 100 |
| V | 2,6-dimethylpyrazine (C6) | 38 | — | — | 56 | — | — | 6 | — | — | 100 |
| V | 2-ethyl-5-methylpyrazine (C7) | 10 | — | — | 22 | 34 | — | — | 45 | — | **111** |
| V | trimethylpyrazine (C7) | 5 | — | — | 12 | 48 | (blank) | 10 | 25 | — | 100 |
| V | furfural (C5) | 29 | — | — | — | 38 | 33 | — | — | — | 100 |
| V | 2-acetylfuran (C6) | 23 | — | — | — | 35 | 15 | 37 | — | — | **110** |
| V | 4-(2-furanyl)-3-buten-2-one (C8) | 12 | — | — | — | — | 21 | — | 13 | 54 | 100 |
| VI | 2-methylthiophene (C5) | 42 | — | — | — | — | 48 | — | — | — | **90** |
| VI | 2-methyl-3-furanthiol (C5) | 25 | — | — | — | 34 | 41 | — | — | — | 100 |
| VI | ethylpyrazine (C6) | 15 | (blank) | (blank) | 43 | (blank) | 26 | 23 | — | — | **107** |
| VI | 2-furfurylthiol (C5) | 52 | — | — | — | — | 48 | — | — | — | 100 |
| VII | 5-ethenyl-4-methylthiazole (C6) | 100 | — | — | — | — | — | — | — | — | 100 |
| VII | 2-furanmethanol (C5) | 48 | — | — | — | 20 | 32 | — | — | — | 100 |

Group V peptide + xylose/[13C5]-xylose; VI + Cys; VII + thiamine. Group VIII (all four) has no rows.

**Derived [Z]** (expected under a 1:1 pool with no isotope effect: one intact C5 → M 50 / M+5 50; two independent
fragments → M 25 %):

| compound | unlabelled share as printed | if M+4 = [M−H]+ of M+5, plus matching unseen [M−H]+ of M | reading |
|---|---|---|---|
| furfural | 29 % | 47 % (r = 38/33) | intact C5, by construction (xylose is its only source) |
| 2-methylfuran | 36 % | 44 % (r = 19/45) | intact C5 |
| MFT | 25 % | 38 % (r = 34/41) | intact C5; no excess unlabelled carbon |
| 2-furanmethanol | 48 % | 60 % (r = 20/32) | intact C5 (or partly unlabelled source; ±10 points noise) |
| FFT | 52 % | — (no M+4) | intact C5 |
| 2-methylthiophene | 47 % of 90 | — | intact C5 |
| 5-ethenyl-4-methylthiazole | 100 % | — | no xylose carbon: thiamine-derived |

The spread of the unlabelled share among compounds that must be intact-C5 (38-60 %) sets the noise floor of these
numbers at about ±10 percentage points [Z]. 2-acetylfuran M+6 37 % and 4-(2-furanyl)-3-buten-2-one M+8 54 % require
an extra labelled C1 / C3 at rates well above the 25 % a 1:1 pool allows for a second independent unit [Z]; with no
spectra shown, these rows cannot be checked.

### 4c. Fig. 2 (p. 3): furfural, Group V labelled sample, ×10² counts, read ±0.02 [F]

| m/z | 95 | 96 | 97 | 99 | 100 | 101 |
|---|---|---|---|---|---|---|
| intensity | 0.53 | 0.55 | ~0.04 | ~0.08 | 0.98 | 0.93 |

NIST reference (upper panel): 95 ≈ 96 ≈ 100 % [F]. [Z]: (95+96) : (99+100+101) = 1.08 : 1.99 = **35 : 65**; the
authors' m/z 100 "[M+ +4]" equals the [M−H]+ of the M+5 isotopologue. Table 3's 29 / 38 / 33 does not reproduce from
these apex intensities (they give 22 / 40 / 38 dropping m/z 95), so Table 3 is probably from another scan or from ion
chromatograms; basis unstated.

### 4d. Fig. 1 (p. 3): ethylpyrazine (MW 108), Group VI labelled sample, read ±0.03 axis units [F]

| m/z | 107 | 108 | 109 | 110 | 111 | 112 | 113 | 114 |
|---|---|---|---|---|---|---|---|---|
| intensity | 0.51 | 0.45 | 0.60 | 0.50 | 1.26 | 0.66 | 1.20 | 0.20 |

NIST reference: 107 ([M−H]+) 100, 108 (M+•) ~72 [F]. The authors call 107 "[M+]" and 111 / 112 / 113 "M+4 / M+5 /
M+6" (Fig. 1) but tabulate M+3 / M+5 / M+6 (Table 3). Deconvolving with a constant [M−H]+ : M+• = 1 : 0.72 gives a
negative M+5 abundance [Z]; the cluster is not one compound with one fragmentation ratio (co-elution with
2,6-dimethylpyrazine, same nominal mass 108, RI 1306 vs 1314, P ~0.5). **Ethylpyrazine row unusable.**

### 4e. Fig. 3 (p. 4): 5-ethenyl-4-methylthiazole (MW 125), Group VII, ×10⁻¹, read ±0.05 [F]

m/z 124 ≈ 0.75, 125 ≈ 4.1, 126 ≈ 0.70, 127 ≈ 0.2. 126/125 ≈ 17 % vs ~8 % expected from natural 13C/33S/15N [Z];
a small labelled contribution cannot be excluded (P ~0.3), but no M+3..M+6 cluster is visible: the authors' "not
detected in the labeled xylose carbon skeleton" holds within read error.

## 5. Flags and internal consistency

1. [!] **[M−H]+ counted as an isotopologue** (Fig. 2 furfural "M+4" = 100; Fig. 1 ethylpyrazine "M+" = 107). No
   natural-abundance correction; no M+1 reported anywhere.
2. [!] **Row sums 111 (2-ethyl-5-methylpyrazine), 110 (2-acetylfuran), 107 (ethylpyrazine), 90
   (2-methylthiophene).**
3. [!] **Text ratios contradict Table 3** (methylpyrazine "1:1:1:3" vs 18:23:25:34; 2,6-DMP "1:2:1" vs 38:56:6;
   2-ethyl-5-methylpyrazine and trimethylpyrazine "1:2:2:1" vs 10:22:34:45 and 5:12:48:(M+6 10):25). Text says
   furfurylthiol and 2-furanmethanol show "M+, M+5 and M+6"; the table has no M+6 for either and an M+4 for
   2-furanmethanol. Text calls furfural's M+4 and M+5 both "singly labeled ([M++4])".
4. [!] **Thiazole identity**: abstract, Fig. 3 caption and Fig. 6 say 4-methyl-5-hydroxyethylthiazole (MW 143); the
   spectrum (m/z 125), Table 2 and Table 3 are 5-ethenyl-4-methylthiazole.
5. [!] **Labelled pots have twice the sugar** of their unlabelled twins (100 vs 50 µmol); furfural and 2-acetylfuran
   carry Table 3 isotopologues in Group V but are absent ("—") from Group I in Table 2.
6. [!] **Table 2 anomalies**: 2,5-dimethylfuran 20,742,697.3 in Group I (misprint or wrong peak); hexanal 1380 and
   octanal 45 in Group IV, which contains no lipid; Group IV numbers 10-1000 × the others with no explanation;
   2-ethyl-3,5- and 3-ethyl-2,5-dimethylpyrazine share RI 1427 (likely one peak); RI order broken (1314 before 1309).
7. [!] **Counts**: text says 31 compounds in the YE pot and 12 in Group I; Table 2 lists 29 and 10 [Z]; "8 furans"
   vs 10 furan rows.
8. [!] **IS concentration** stated two ways (500 ng/g vs ~102 ng/mL); no unit on Table 2.
9. [!] **Citation errors**: "Wang, Yang & Song 2017, Food Chem 131:280" is Wang 2012 (`wang2012_extraction.md`);
   "Christoph & Matthieu 2007" and "Christoph & Tomas 2003" are Cerny & Briffod 2007 and Cerny & Davidek 2003 (given
   names as surnames); "Monika & Schieberle 2009" duplicates Christlbauer & Schieberle 2009.
10. [!] **Title vs content**: "from yeast extract using CAMOLA", but CAMOLA was run only on synthetic-peptide pots.

## 6. Existing records in the repo that cite this paper

Grep 2026-11-03 (`10.1016/j.lwt.2019.06.006`, `alim2019`, `Alim 2019`, `Alim et al`, `CAMOLA` scanned for Alim) over
`data/`, `docs/`, `src/`, `scripts/`, `tasks/`, `results/`: **no ledger, registry, fit row or code record**. Three
text mentions:

| file | field / claim | paper | verdict |
|---|---|---|---|
| `docs/guides/CITATION_CONTEXT_2026-10.md` l. 135 | title, DOI | same | **match** |
| same | class B, "yeast-extract meat-flavour background; retention-index method" (how citers use it) | uses LRI on DB-WAX; YE background | **match** as a description of citers |
| same | representative context "...mainly derived from sulfur-containing nitrogen compounds" | a citer's sentence; not in the paper | **not a paper claim** (correctly presented as a citer context) |
| same | recommendation SKIP | paper does hold label tables, but small, uncorrected, one condition | **defensible**; the label data are a weak check only |
| `docs/guides/LITERATURE_SEARCH_2026-10.md` l. 151 | "CAMOLA labelling of beef-meaty odorants from yeast extract" | labelling was in synthetic-peptide + xylose pots, never in yeast extract | **mismatch** (the YE pot was unlabelled) |
| `data/keys/pdf_dossier_aliases.yml` l. 93 (triage, not edited) | "from yeast extract; route shares (E6, E13)" | no YE labelling; no MFT thiamine-vs-sugar share (Group VIII unreported); pyrazine C-module patterns only | **partly mismatched**: gives intact-C5 evidence without thiamine, not route shares |

## 7. Proposed ledger entries (NOT written into the ledgers)

```yaml
# sulfur.yml
  - id: alim2019_table3_group_vi_camola_xylose_isotopomers
    dossier: alim2019_extraction.md
    paper: "Alim, Song, Liu, Zou, Zhang, Zhang & Raza 2019 LWT 112:108239"
    lab: "BTBU (Song H)"
    system: "synthetic YE peptide (Pro-Glu-Thr / Glu-Asn / Arg-Leu, which not stated) 100 umol + xylose 50 + [13C5]-xylose 50 umol + L-cysteine 50 umol, no thiamine (Group VI); volume not stated"
    medium: "water, 'pH 6.0' (phosphate implied, strength not stated), vial sealed with aluminium foil inside a Parr high-pressure reactor"
    temperature_c: [135]
    ph: {initial: 6.0, measured_during: false, final: null}
    time_min: [90]
    species: ["MFT", "FFT", "2-methylthiophene", "ethylpyrazine"]
    kind: ratio
    location: "Table 3 (isotopologue %, Group VI rows); Fig. 1 ethylpyrazine spectrum"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "the following conditions: 90 min, 135 °C and pH 6.0."
    flags: [condition_uncertain, ph_start_only, volume_not_stated, no_mh_correction, no_natural_abundance_correction, "row_sum_90: 2-methylthiophene", "row_sum_107: ethylpyrazine", "complements: alim2019_table3_group_v_vii_camola_pyrazines_furans"]
    note: "MFT 25/34/41 (M/M+4/M+5), FFT 52/48, 2-methylthiophene 42/48: intact xylose C5 in a no-thiamine pentose + Cys pot (agrees with cerny2007b system F). M+4 is probably [M-H]+ of M+5 (shown for furfural, Fig. 2); MFT re-allocated 38:62. No thiamine share: labelled Group VIII has no rows. Ethylpyrazine row on the [M-H]+ mass axis: do not use. Labelled pot has 2x the xylose of unlabelled Group II. BTBU Song group, not independent of wang2012 / xu2013."
  - id: alim2019_table2_sulfur_rows
    dossier: alim2019_extraction.md
    paper: "Alim, Song, Liu, Zou, Zhang, Zhang & Raza 2019 LWT 112:108239"
    lab: "BTBU (Song H)"
    system: "Groups I-IV: peptide 100 + xylose 50 umol +/- Cys 50 +/- thiamine 50 umol; and a YE paste LA00 + xylose + Cys + thiamine pot (amounts not stated)"
    medium: "water, 'pH 6.0', foil-sealed vial in Parr reactor; YE column is a food matrix"
    temperature_c: [135]
    ph: {initial: 6.0, measured_during: false, final: null}
    time_min: [90]
    species: ["MFT", "FFT", "2-methylthiophene", "5-ethenyl-4-methylthiazole"]
    kind: endpoint
    location: "Table 2 (sulfur rows)"
    status: digitised_in_dossier
    quant_class: semi_quantitative
    simulable: false
    condition_quote: "The test was repeated three times under the same conditions."
    flags: [condition_uncertain, ph_start_only, volume_not_stated, unit_not_printed, no_sd, single_is_no_response_factor, matrix_food, "is_concentration_inconsistent: 500 ng/g vs 0.816 ug in 8 mL", "complements: alim2019_table2_pyrazine_furan_rows"]
    note: "IS ratio (2-methyl-3-heptanone, RF = 1), unit not printed. MFT / FFT / 2-MT only with Cys (II, IV); thiamine raises MFT 5.83 -> 18.45 and FFT 7.71 -> 249.12 (II -> IV), directional only. YE recipe and composition not printed (maybe Supplementary Table 1)."

# aroma_other.yml
  - id: alim2019_table3_group_v_vii_camola_pyrazines_furans
    dossier: alim2019_extraction.md
    paper: "Alim, Song, Liu, Zou, Zhang, Zhang & Raza 2019 LWT 112:108239"
    lab: "BTBU (Song H)"
    system: "synthetic YE peptide 100 umol + xylose 50 + [13C5]-xylose 50 umol (Group V); + thiamine 50 umol (Group VII); volume not stated"
    medium: "water, 'pH 6.0', foil-sealed vial in Parr reactor"
    temperature_c: [135]
    ph: {initial: 6.0, measured_during: false, final: null}
    time_min: [90]
    species: [methylpyrazine, "2,5-dimethylpyrazine", "2,6-dimethylpyrazine", "2-ethyl-5-methylpyrazine", trimethylpyrazine, furfural, 2-methylfuran, 2-acetylfuran, "4-(2-furanyl)-3-buten-2-one", 2-furanmethanol, "5-ethenyl-4-methylthiazole"]
    kind: ratio
    location: "Table 3 (Group V and VII rows); Fig. 2 furfural spectrum (digitised, dossier sec. 4c)"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "A 1:1 fully labeled and unlabeled D-xylose, was selected"
    flags: [condition_uncertain, ph_start_only, volume_not_stated, no_mh_correction, no_natural_abundance_correction, text_ratios_contradict_table, "row_sum_111: 2-ethyl-5-methylpyrazine", "row_sum_110: 2-acetylfuran", "complements: alim2019_table3_group_vi_camola_xylose_isotopomers"]
    note: "Pyrazine carbon from xylose C2/C3 modules: methylpyrazine 18/23/25/34 (M/M+2/M+3/M+5) and 2,5-DMP 22/48/30 (M/M+3/M+6) near binomial 1:1 expectations; 2,6-DMP, 2-ethyl-5-methylpyrazine, trimethylpyrazine fit no simple module model. No N tracing. Furan M+4 = [M-H]+ of M+5 (Fig. 2). Thiazole 100 % unlabelled (thiamine). Labelled pots carry 2x the xylose of Groups I-IV."
  - id: alim2019_table2_pyrazine_furan_rows
    dossier: alim2019_extraction.md
    paper: "Alim, Song, Liu, Zou, Zhang, Zhang & Raza 2019 LWT 112:108239"
    lab: "BTBU (Song H)"
    system: "Groups I-IV (peptide + xylose +/- Cys +/- thiamine) and a YE paste pot, as alim2019_table2_sulfur_rows"
    medium: "water, 'pH 6.0', foil-sealed vial in Parr reactor; YE column is a food matrix"
    temperature_c: [135]
    ph: {initial: 6.0, measured_during: false, final: null}
    time_min: [90]
    species: [pyrazines (9), furans, 2-methylbutanal, benzaldehyde, hexanal, octanal, nonanal, 2-acetylpyrrole, indole]
    kind: endpoint
    location: "Table 2 (non-sulfur rows)"
    status: digitised_in_dossier
    quant_class: semi_quantitative
    simulable: false
    condition_quote: "The test was repeated three times under the same conditions."
    flags: [condition_uncertain, ph_start_only, volume_not_stated, unit_not_printed, no_sd, single_is_no_response_factor, matrix_food, "cell_outlier: 2,5-dimethylfuran Group I 20742697.3", "lipid_aldehydes_in_lipid_free_pot: Group IV hexanal 1380", "complements: alim2019_table2_sulfur_rows"]
    note: "Directional only. 3-ethyl-2,5-DMP dominates every peptide pot (818-33465); 2-ethyl-3,5- and 3-ethyl-2,5-DMP share RI 1427. Text counts (31 YE compounds, 12 in Group I) do not match the table (29, 10)."
```

Why no fit-grade entry: no concentrations (volume unstated), peptide reactants not in the engine, one T / t, no SD,
uncorrected ion basis. The best-value cell is the Group VI MFT / FFT intact-C5 result as a directional check
(P ~0.3 it is worth wiring as a hold-out, given Cerny 2007b system F already covers the same claim with a cleaner
design).

## 8. Supplementary material

Yes: "Appendix A. Supplementary data" online at the DOI, and sec. 2.6 says the tests were designed "as shown in
Supplementary Material Table 1". **Not on disk** (`data/articles/` has only `alim2019.pdf`). It may hold the YE-pot
recipe (absent from the paper) and possibly volumes; it is unlikely to hold spectra or corrected isotopologue data
(P ~0.2). Fetching it is low priority: it would not make any entry fit-grade unless it prints the water volume and
the Group VIII labelled data (P ~0.1).
