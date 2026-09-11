# Huang, Tippmann & Becker — EXTRACTION (glucose or maltose + leucine or isoleucine, 10:1 molar, 0.1 mol/L phosphate pH 5.2, 90-130 C, 0-360 min; the SAME paper already extracted as `huang2017_extraction.md`)

### `data/articles/huang2016.pdf` and the source of `huang2017_extraction.md` are one file and one paper; there is no second Huang paper on disk, so nothing here fixes the missing concentrations — the pot's absolute charges are still not printed, and the constants still do not convert.

**Source on disk:** `data/articles/huang2016.pdf` (11 pp., owner's download, 2026-09-08; md5
`d4a04696afe4786d36bdfe04fb9fcea8`). It is the **only** file matching `huang*` in `data/articles/`.
Read from a fresh `pdftotext -layout` extraction into the scratchpad and searched again, whole-file, for
any absolute concentration: there is none. Tables 1 and 2 came through the text layer with the exponent
signs lost (`5.68E203` for 5.68E-03, `e22953.8/T` for e^(-2953.8/T)); each sign was re-derived from the
printed activation energy below and is not in doubt. Figures 1-7 are FIGURE-ONLY and their axis units
are not in the text; no value was read off any of them.

## 0. Identity

| field | value |
|---|---|
| Title | "A kinetic study on the formation of 2- and 3-methylbutanal" |
| Authors | Yarong Huang, Johannes Tippmann, Thomas Becker (Institute of Brewing and Beverage Technology, TU München, Weihenstephaner Steig 20, 85354 Freising) |
| Venue | Journal of Food Process Engineering (Wiley). Received 14 Aug 2015, accepted **3 Feb 2016**; the PDF's running footer still reads "00 (2016) 00-00" (an accepted-article proof), and the article issued in the 2017 volume as e12375 |
| DOI | 10.1111/jfpe.12375 |
| Funding | FEI / AiF project 16968N |
| **Duplicate status** | **This is the paper already extracted in `data/lit/extraction_dossiers/huang2017_extraction.md`.** That dossier's own source line names this same file: "`data/articles/huang2016.pdf` ... the file name says 2016 because the article was accepted February 2016 and issued in J. Food Process Eng. 2017". Title, authors, DOI, both tables and every number match. There is no separate 2016 paper by this group on disk (flag 1) |
| Naming | 2-MB = 2-methylbutanal (from isoleucine); 3-MB = 3-methylbutanal (from leucine); GL / ML / GI / MI = glucose-leucine / maltose-leucine / glucose-isoleucine / maltose-isoleucine; ZP = the lumped unmeasured "intermediate" pool; AP = the aldehyde's "degradation products" |

## 1. Why it matters

`results/validation/kinetic_core_b19_prereg_draft.md` section 5 refuses a per-amino-acid Strecker rate
wave because no per-amino-acid rate exists in the corpus, and leaves an identity-ratio layer on
glycine's fitted step (`FROZEN_B18` in `src/kinetic_core/parameters_pyrazine.py`) as the only honest
form. `huang2017_extraction.md` was read for that wave and found unusable for one reason, its flag 1:
the paper prints only "mol ratio 10:1" and "0.1 mol/L phosphate", never an absolute sugar or amino-acid
concentration, so its pseudo-first-order `k1` cannot be turned into the second-order (dicarbonyl x amino
acid) form that `FROZEN_B18` uses, and its aldehyde figures cannot be converted to mmol/L.

This dossier was commissioned to check whether a 2016 paper by the same group repairs that. **It does
not, and the reason is not a defect in the 2016 paper: there is no 2016 paper.** The file called
`huang2016.pdf` is the accepted-article proof of the same JFPE article, and its two tables are the two
tables already transcribed in `huang2017_extraction.md`. Re-reading it whole confirms the refusal
stands: the only quantitative statements about how much was charged are "**mol ratio 10:1**" and "**a
phosphate buffer (0.1 mol/L; pH 5.2)**", plus a justification for the ratio ("the molar ratio between
sugar and amino compound of wort is more than 3:1 (Fox et al. 1983)"). No mass, no molarity, no volume,
no headspace.

**Read `huang2017_extraction.md` as the primary record for this paper.** This file exists so that a
future reader who opens `huang2016.pdf` is not misled into thinking it is a second source, and so that
the question "does the 2016 paper print the pot's concentrations, so its constants convert?" has a
recorded answer: **no**.

## 2. Methods as they matter to a model

Unchanged from `huang2017_extraction.md` section 2 and re-verified against the file. In brief, and only
the parts that bear on the concentration question:

- **Charge (verbatim, the whole of it):** "Disaccharide sugar (maltose) or monosaccharide sugar
  (glucose) and amino acid (leucine or isoleucine), **mol ratio 10:1**, were dissolved in a **phosphate
  buffer (0.1 mol/L; pH 5.2)**". Every subsequent mention of amount in the paper is a ratio or a
  percentage of the initial value. The whole-file search for `mol/L`, `mmol`, `g/L`, `mg/mL`, `g of`
  returns only: the 0.1 mol/L buffer; the 400 mmol/L NaOH of the sugar HPLC eluent; the 1:10 dilution
  of the amino-acid samples; the 100 µL + 900 µL sugar injection mix; and the GC's 4 µL injection. None
  of those is a reactant concentration.
- **Buffer / pH:** 0.1 mol/L phosphate, pH 5.2 (wort pH). Whether the pH was followed during heating is
  not stated.
- **Vessel and heating:** "closed glass tubes" in a "heating block", "various times (**0-360 min**) at
  different temperatures (**90-130 °C**)", cooled in ice water at the planned time. From the Results
  the ladder is 90, 100, 110, 120, 130 C. **Tube volume and headspace are not stated**, which matters
  because 2-MB and 3-MB boil at 91 C (the paper says so itself).
- **Replicates:** "The tests were repeated three times." The only error figures printed are the ± on
  the Table 1 constants, and the paper does not say whether they are a standard deviation, a standard
  error or a regression interval.
- **Analytics:** sugars by HPAEC-PAD (Dionex ICS-1000, CarboPac PA10; internal standard not named);
  amino acids by OPA/FMOC derivatisation and HPLC (Dionex UltiMate 3000); 2-MB and 3-MB by
  purge-and-trap GC-FID (HP 5890, Chrompack PTI, HP-Innowax and HP Ultra-2, purge vessel 50 C, hydrogen
  carrier, FID 250 C), "following the method of MEBAK with major modification (Pfenninger 1993)".
  **No calibration, internal standard, LOD or recovery is given for the aldehydes**, and the
  concentration unit of the aldehyde figures is not in the text.
- **Model (Scheme 3, Eqs. 9-12):** d[Leu]/dt = -k1[Leu]; d[ZP]/dt = k1[Leu] - k2[ZP];
  d[3MB]/dt = k2[ZP] - k3[3MB]; d[AP]/dt = k3[3MB]. First order in each species; **the sugar does not
  appear in the rate law at all** (it is in ten-fold excess and its loss is described separately as
  linear in time). Fitted by multiresponse regression in Athena Visual Studio. k3 "is the smallest and
  has the largest corresponding interval" and is not tabulated.
- **Statistics:** Scheffé F-test between sugar/amino-acid pairs, p < 0.05.

## 3. Tables re-typed

The paper contains exactly two tables. Both are re-typed here as printed, so this file stands alone;
they are identical to the two in `huang2017_extraction.md`.

### Table 1. "The reaction rate constants of the reaction between GL, ML, GI and MI at 100 °C"

Units as printed: "/min" (min^-1). Values are mean ± (basis of the ± not stated).

| constant | GL | ML | GI | MI |
|---|---|---|---|---|
| k1 (/min) | 5.68E-03 ± 6.74E-04 | 7.88E-04 ± 1.39E-05 | 1.99E-04 ± 2.03E-05 | 3.24E-04 ± 9.82E-05 |
| k2 (/min) | 1.14E-03 ± 9.09E-05 | 1.62E-03 ± 6.76E-04 | 3.61E-04 ± 7.90E-05 | 1.20E-03 ± 2.20E-04 |

**Amino-acid half-lives at 100 C from k1 (mine):** 122 min (GL), 880 min (ML), 3483 min (GI), 2139 min
(MI). **In s^-1 (mine):** k1 = 9.47e-5, 1.31e-5, 3.32e-6, 5.40e-6; k2 = 1.90e-5, 2.70e-5, 6.02e-6,
2.00e-5.

### Table 2. "Kinetic parameters with temperature dependence for the model in Schema 3"

Re-typed as k = A·exp(-B/T), T in kelvin, k in min^-1 (the signs are restored from the printed Ea; the
text layer prints them as `e22953.8/T` etc.).

| row | GL | ML | GI | MI |
|---|---|---|---|---|
| Arrhenius-function 1 (k1) | k = 3.015 · e^(-2953.8/T) | k = 2.036 × 10³ · e^(-5448.8/T) | k = 3.895 × 10⁴ · e^(-7000.8/T) | k = 1.622 × 10³ · e^(-5732.6/T) |
| Arrhenius-function 2 (k2) | k = 1.434 × 10⁹ · e^(-10020/T) | k = 6.392 × 10¹² · e^(-13283/T) | k = 1.892 × 10¹⁴ · e^(-14515/T) | k = 2.651 × 10¹² · e^(-13259/T) |
| R²₁ | 0.9891 | 0.9514 | 0.9355 | 0.9816 |
| R²₂ | 0.9874 | 0.9664 | 0.9369 | 0.9704 |
| Ea1 (kJ/mol) | 24.56 | 45.3 | 58.2 | 47.67 |
| Ea2 (kJ/mol) | 83.31 | 110.43 | 120.68 | 110.24 |
| Etotal (kJ/mol) | 107.87 | 155.73 | 178.88 | 157.91 |

**Arithmetic re-checked independently for this dossier (mine).** Slope × R: 2953.8 × 8.314 = 24.56;
5448.8 × 8.314 = 45.30; 7000.8 × 8.314 = 58.20; 5732.6 × 8.314 = 47.66; 10020 × 8.314 = 83.31;
13283 × 8.314 = 110.43; 14515 × 8.314 = 120.68; 13259 × 8.314 = 110.24 kJ/mol. Every printed Ea is the
slope of its own printed line, which fixes the lost signs. "Etotal" is exactly Ea1 + Ea2 and is not the
barrier of any measurable quantity.

**The two tables disagree for the glucose systems (mine, reproducing `huang2017_extraction.md` flag 2).**
Evaluating each line at 373.15 K and dividing the Table 1 value by it:

| constant | GL | ML | GI | MI |
|---|---|---|---|---|
| k1 from line at 100 C (min^-1) | 1.10e-3 | 9.27e-4 | 2.77e-4 | 3.45e-4 |
| Table 1 k1 ÷ line | **5.16** | 0.85 | 0.72 | 0.94 |
| k2 from line at 100 C (min^-1) | 3.12e-3 | 2.22e-3 | 2.42e-3 | 9.81e-4 |
| Table 1 k2 ÷ line | **0.36** | 0.73 | **0.15** | 1.22 |

Lines evaluated over the ladder (min^-1), for the record: GL k1 8.85e-4 (90 C) → 1.98e-3 (130 C);
GL k2 1.49e-3 → 2.30e-2; ML k1 6.20e-4 → 2.75e-3; ML k2 8.32e-4 → 3.14e-2; GI k1 1.65e-4 → 1.12e-3;
GI k2 8.29e-4 → 4.37e-2; MI k1 2.26e-4 → 1.08e-3; MI k2 3.69e-4 → 1.38e-2.

### Statements in the text (the curves are FIGURE-ONLY)

- Leucine with glucose: >60 % degraded after 360 min at 130 C, about 20 % at 90 C.
- Isoleucine with maltose: 35 % after 360 min at 130 C, about 10 % at 90 C, "constant at the beginning
  or even increased slightly at 90 °C"; with glucose at 110-130 C the loss "seem[s] to occur in three
  steps".
- "Compared to isoleucine, leucine is degraded about **two times faster** than isoleucine."
- Glucose: about 40 % lost after 360 min at 130 C with leucine, about 30 % with isoleucine; no
  significant difference at 90-120 C. Mannose, fructose, maltose, maltotriose and sucrose were detected
  in the glucose samples. Maltose loss is linear and "a little bit slower" than glucose's.
- 3-MB: raising the temperature from 90 to 130 C increases the concentration "**about 100 times**"; the
  90 and 100 C curves are close, the 100 C values significantly higher.
- Ordering: "**GL > ML > GI > MI**"; "the type of sugar ... influences the amount of Strecker flavors
  greater than the type of amino acids."
- Secondary quotations (other people's numbers, not measured here): Cremer & Eichner 2000, zero order
  over 10-120 min, Ea 124 (2-MB) and 120 (3-MB) kJ/mol; Chan & Reineccius 1994, pseudo-zero order,
  Ea(3-MB) 80.4 kJ/mol; Balagiannis et al. 2009, first two steps 137 ± 15.2 and 48.7 ± 8.4 kJ/mol;
  Gomyo et al. 1989, blue pigment Ea 73 kJ/mol.

## 4. Kinetic numbers the repository can use

Registry (`data/keys/compounds.yml`): `3_methylbutanal` and `2_methylbutanal` are present. Leucine,
isoleucine, glucose, maltose and the lumped intermediate have no id (`data/lit/reaction_rules.yml` uses
the short names Leu, Ile, Glc).

**Every row below is already carried in `huang2017_extraction.md` section 4.** They are repeated here
only so this file is self-contained; they are one set of numbers, not two, and must not be fitted twice.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| k1, amino-acid loss, GL / ML / GI / MI | 5.68e-3 ± 6.74e-4 / 7.88e-4 ± 1.39e-5 / 1.99e-4 ± 2.03e-5 / 3.24e-4 ± 9.82e-5 | min^-1 | 0.1 mol/L phosphate pH 5.2, sugar : amino acid 10:1, **absolute concentrations not printed**, 100 C, closed tubes, 0-360 min, n = 3 | Table 1 | measured_rate (lumped loss, all sinks; pseudo-first order) |
| k2, intermediate → aldehyde, GL / ML / GI / MI | 1.14e-3 ± 9.09e-5 / 1.62e-3 ± 6.76e-4 / 3.61e-4 ± 7.90e-5 / 1.20e-3 ± 2.20e-4 | min^-1 | same, 100 C | Table 1 | measured_rate (model-identified only; the intermediate was never measured) |
| Ea1, GL / ML / GI / MI | 24.56 / 45.3 / 58.2 / 47.67 | kJ/mol | 90-130 C, five temperatures, R² 0.94-0.99 | Table 2 | measured_barrier (of the lumped loss) |
| Ea2, GL / ML / GI / MI | 83.31 / 110.43 / 120.68 / 110.24 | kJ/mol | same | Table 2 | measured_barrier (see the two-table disagreement above) |
| Arrhenius prefactors A1 / A2 | 3.015, 2.036e3, 3.895e4, 1.622e3 / 1.434e9, 6.392e12, 1.892e14, 2.651e12 | min^-1 | same | Table 2 | measured, printed |
| k1(GL)/k1(GI) and k1(ML)/k1(MI) at 100 C | 28.5; 2.43 | – | same | Table 1 (mine) | within_study_ratio (leucine vs isoleucine; the text says "about two times", which matches only the maltose pair) |
| k1(GL)/k1(ML) and k1(GI)/k1(MI) | 7.2; 0.61 | – | same | Table 1 (mine) | within_study_ratio (glucose vs maltose; the sign reverses between the two amino acids) |
| aldehyde formation ordering | GL > ML > GI > MI | – | all temperatures | text | within_study_ratio (ordinal only) |
| leucine / isoleucine / glucose losses at 360 min | >60 % and ~20 % (Leu, GL, 130 / 90 C); 35 % and ~10 % (Ile, MI); ~40 % and ~30 % (glucose, GL / GI at 130 C) | % | as stated | text | level_only (approximate) |
| 3-MB at 130 C vs 90 C | "about 100 times" | – | GL | text | within_study_ratio (the authors' reading of their own figure) |
| "Etotal" 107.87-178.88 | – | – | – | Table 2 | **not a barrier**: an arithmetic sum of two step barriers; do not carry |
| k3, aldehyde loss | not tabulated | – | – | text | not reported |
| any absolute concentration of sugar, amino acid, intermediate or aldehyde | **none anywhere in the paper** | – | – | whole file searched | not reported |
| 2-MB, 3-MB, amino acid, sugar vs time | – | unit not in the text | 90-130 C | Figs. 1-7 | figure_only |

### Does the 2016 paper fix what `huang2017_extraction.md` could not use?

**No, and it cannot, because it is the same paper.** `huang2016.pdf` is the accepted-article proof of
JFPE e12375, the document `huang2017_extraction.md` was written from; the dossier's own source line says
so. A second reading of the whole file adds no charge, no molarity, no mass and no volume. Therefore:

- **k1 still does not convert.** It is a pseudo-first-order loss constant of the amino acid at an
  unknown concentration with an unknown ten-fold sugar excess. Turning it into the second-order form
  `FROZEN_B18` uses — a constant on [dicarbonyl] × [amino acid], in L/(mmol·min) — needs the sugar
  concentration, which is not printed. Neither can a molar yield be recovered, because the aldehyde
  figures carry no unit.
- **k1 is not a Strecker rate in any case.** It lumps Amadori formation, the Strecker step and
  melanoidin binding, as `huang2017_extraction.md` flag 3 already records; the 24.6 kJ/mol barrier of
  the GL system is what a diffusion- or equilibrium-limited lumped step looks like, not a bond-breaking
  barrier.
- **What survives is one within-study identity ratio**, leucine against isoleucine at fixed sugar and
  fixed pot: 28.5 on glucose and 2.43 on maltose from Table 1, against the text's "about two times".
  A ratio that moves by an order of magnitude between two sugars in the same paper is a weak row, and
  `huang2017_extraction.md` flag 4 already says to carry both values or neither.
- **What a fix would take:** the pot's absolute concentrations, from the authors or from the underlying
  thesis. That request is flag 3 below.

## 5. Flags

1. **Duplicate source.** `data/articles/huang2016.pdf` is the paper already extracted as
   `huang2017_extraction.md`; there is no second Huang paper in `data/articles/`. Do not treat this
   dossier and that one as two sources, and do not fit Tables 1 and 2 twice. `huang2017_extraction.md`
   is the primary record; this file records the collision and the answer to the concentration question.
   If the repository wants one file per PDF, this is that file, and it should carry a pointer, not a
   second vote.
2. **No absolute concentrations, confirmed on a second whole-file reading.** Only "mol ratio 10:1" and
   "0.1 mol/L phosphate, pH 5.2". Nothing in this paper converts to mmol/L, and nothing in it can be
   written as a second-order constant.
3. **To request from the authors:** the sugar and amino-acid concentrations of the model solutions, the
   glass-tube volume and headspace, and the concentration unit of Figures 5 and 6. The corresponding
   author's address is printed (huangyarong@hotmail.de); the thesis-level detail may sit with the
   FEI/AiF project 16968N report.
4. **Tables 1 and 2 disagree for the glucose systems** (GL k1 is 5.2x its own Arrhenius line, GI k2 is
   6.7x below its own line, while ML and MI agree within 30 %). The paper never says how the two tables
   were obtained. Neither should be used alone.
5. **"Etotal" is not a barrier.** The abstract's "total activation energies ... ranging from 107.87 to
   178.88 kJ/mol" is Ea1 + Ea2, an arithmetic sum, and must not be carried as an activation energy of
   methylbutanal formation.
6. **k2 exists only through the model.** The intermediate pool ZP was never measured, so k2 and Ea2
   depend entirely on the assumed three-step chain.
7. **Aldehyde quantification undocumented**, and 2-MB / 3-MB boil at 91 C while the ladder runs to
   130 C. The tubes were closed, but the headspace is not stated, so the liquid-phase aldehyde depends
   on an unknown gas/liquid ratio at every temperature above 91 C.
8. **The rate law has no sugar term and no dicarbonyl**, so the constants belong to the pH 5.2, 10:1
   regime only and say nothing about the dicarbonyl supply that `PYRAZINE_SUPPLY_CAVEAT` names as the
   trunk's real uncertainty.
9. **Registry gaps:** leucine, isoleucine, glucose, maltose and the intermediate pool have no
   `compounds.yml` row; the two aldehydes do.
