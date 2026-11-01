# Van Lancker, Adams & De Kimpe 2010 — EXTRACTION
### *"Formation of Pyrazines in Maillard Model Systems of Lysine-Containing Dipeptides"*
### Eight Lys-X dipeptides (X = Gly, Ala, Leu, Lys, Ser, Glu, Phe, Cys) vs the equimolar free amino acids vs free lysine (0.2 and 0.4 M), each with glucose (0.2 M), methylglyoxal (0.02 M) or glyoxal (0.02 M); unbuffered water, initial pH 8, 130 °C, 2 h, closed Pyrex tubes. Pyrazines by SBSE-GC-MS as **peak areas** with printed THEORETICAL SBSE recoveries; 2,5- and 2,6-DMP co-elute. **The Methods companion of Van Lancker 2012 (whose dossier says "not on disk").**

**Read 2026-10-31 from `data/articles/vanlancker2010.pdf`** (9 pp., JAFC 58:2470-2478; born-digital). Tables 1-2 from
`pdftotext -layout` (column count checked: 17 data columns + recovery in every row used); **Table 3 (printed rotated)
read from a 150 dpi render of p. 2474**; Table 4 and Scheme 1 from the text layer. Derived numbers computed in the
`maillard` env.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Authors | Fien Van Lancker, An Adams, Norbert De Kimpe (Ghent University, Dept. of Organic Chemistry) |
| Venue | *J. Agric. Food Chem.* **2010**, 58 (4), 2470-2478; received 6 Nov 2009, accepted 20 Jan 2010 |
| DOI | **10.1021/jf903898t** (printed; Crossref 2026-10-31: same title, vol. 58, pp. 2470-2478) |

## 1. What this paper answers for the model

| need | what Van Lancker 2010 gives | tag |
|---|---|---|
| pyrazine ratios for R (glucose pots) | PZ, MPZ, 2,5(6)-DMP areas in 17 pots (Table 1). **R(area) = 0.29-3.52**; × 0.84 recovery correction → 0.24-2.97. Lys + Cys: DMP not detected | [M] / [Z] |
| response factors | areas only; printed theoretical SBSE recoveries **PZ 0.3 %, MPZ 0.8 %, DMP 2.0 / 1.6 %** (Gerstel calculator, not measured) [!]; recovery-corrected R = R(area) × 0.3 × 1.8 / 0.8² = **0.84** R(area) [Z]; EI response per mole unknown; **2,5-DMP co-elutes with 2,6-DMP** (2,6-DMP is not an aminoacetone dimer, so the "DMP" here overstates the AKM pool and lowers R) [!] | [M] / [Z] |
| MGO pots / GO pots | MGO: PZ detected only in Lys + Ser (0.05), R = 1.52 there; elsewhere PZ "−". GO: DMP absent; MPZ 0.01-0.20 vs PZ 0.11-6.48, i.e. **MPZ / PZ ≈ 0.01-0.05 (area), ~0.004-0.02 molar** in glyoxal-only pots [Z] | [M] / [Z] |
| 13C / 15N | none measured. Cited: Hwang 1994 (alpha-15N lysine + glucose, 180 °C, pH 8.5): PZ and MPZ took ~40 % of their N from the epsilon-amino group, 2,6-DMP and TMP <= 2 % | [C] |
| condensation k, pH law, Ea | none (one T, one t) | — |
| pyrrole vs pyrazine | pyrroles not reported; "1,4-dimethylpyrrolo(1,2-a)pyrazine" appears (tentative) | [M] |
| peptide vs free amino acid | dipeptides: total pyrazine area 1.3-9.8x the free amino acids with glucose (except Lys-Leu), 1.8-55.5x with MGO; less with glyoxal. Driven by 2,5(6)-DMP and TMP; PZ always LOWER with dipeptides; amino-acid-specific pyrazines (Strecker-aldehyde addition) much lower | [M] / [Q] |
| peptide hydrolysis | Strecker aldehyde areas free / dipeptide ≈ 10 (glucose), 100 (MGO), 100-1000 (GO) (Table 4) [Q]; read as limited hydrolysis | [M] / [Q] |
| mechanism for peptides | imine → 1,5-H shift → 4-hydroxy-2-azadiene → hydrolysis gives the alpha-aminoketone and an alpha-ketoamide, no Strecker aldehyde (Scheme 1) | **[P]** |

**Answers.**

1. **No R above 4 in any pot** (area or recovery-corrected) [Z]. Every rival allows R < 4 (h < 1), so these rows do
   not discriminate P1 / P1x / P1y / P2b by the bound. With co-eluting 2,6-DMP in the denominator they are biased low.
2. **But R barely moves while PZ / DMP moves 100-fold** [Z]. Geometric means: free-amino-acid pots R 0.65 at PZ/DMP
   1.36; dipeptide pots R 1.25 at PZ/DMP 0.013. Under pairing, A/M scales as sqrt(PZ/DMP), so A/M fell ~10x while R
   rose 1.9x. An extra MPZ route whose flux does NOT scale with the aminoacetaldehyde pool A (e = extra / pairing ∝
   1/A) would then need **e ≈ 0.04 in the free-Lys pots** (e ≈ 0.45 in the dipeptide pots) at equal h [Z]. So an
   A-independent extra route contributes ≲ 5-10 % of MPZ in lysine + glucose at 130 °C. An extra route that itself
   scales with A (e.g. A + a hydroxylated C3, P2b's first example) leaves R unchanged and is **not** constrained.
   Confidence in this inference ~0.5 (h may differ between free and peptide pots; 2,6-DMP co-elution; pH drift
   differs; peak areas).
3. **Does it distinguish P1x / P1y from P2b? No.** Point 2 trims only the A-independent flavour of P2b.
4. **Cross-product leak in single-dicarbonyl pots** [Z]: glyoxal-only pots make MPZ at ~1-2 % of PZ (molar), MGO-only
   pots make MPZ at ~0.3-10 % of DMP (area). Any route that turns the aminoacetaldehyde DHP into MPZ by a C1 addition
   (P2a-type) is therefore small in GO + lysine at 130 °C. The authors attribute the MPZ in MGO pots to amino-acid
   carbon (serine) or MGO fragmentation [P].

## 2. Conditions [M]

| item | value |
|---|---|
| Amine | 1 mmol dipeptide, or 1 mmol EACH free amino acid (so the free arm has 2x the alpha-amines), in 5 mL: **0.2 M**; free Lys alone 0.2 and 0.4 M |
| Carbonyl | glucose 1 mmol (**0.2 M**) or MGO / GO 0.1 mmol (**0.02 M**, "to limit self-condensation") |
| Salts | peptides as HBr, HCl or acetate salts (Lys-Ala·HBr, Lys-Lys·HCl, Lys-Leu acetate, ...) [!] counter-ions differ by peptide |
| pH | adjusted to 8 with 2 N NaOH / HCl; **unbuffered** (phosphate and citrate avoided as catalysts); glucose pots fell to **pH 4-5** during heating; MGO / GO pots "remained quite stable" [Q]; re-adjusted to pH 8 on ice before extraction |
| T / t / vessel | 130 °C oil bath, 120 min, 10 mL Pyrex tubes "tightly closed"; ice quench. Headspace / O2 not controlled |
| Extraction | 4.5 mL in a 15 mL vial, 30 min at 35 °C equilibration, SBSE (PDMS Twister 10 mm x 0.5 mm) 30 min at 35 °C; TDS-CIS; HP5-MS; areas x 1e8. Repeatability on standards: SD < 11 % (triplicate). No replicates of the pots stated [!] |

## 3. Quantitative results

### 3.1 Table 1 (glucose pots), the three species used for R [M]; R and ratios [Z]

| pot | PZ | MPZ | 2,5(6)-DMP | R (area) | R x 0.84 | PZ / DMP |
|---|---:|---:|---:|---:|---:|---:|
| Lys 0.2 M | 3.04 | 1.82 | 0.85 | 1.28 | 1.08 | 3.58 |
| Lys 0.4 M | 5.30 | 3.06 | 2.84 | 0.62 | 0.52 | 1.87 |
| LysLys | 0.47 | 1.26 | 11.81 | 0.29 | 0.24 | 0.040 |
| Lys + Gly | 1.63 | 2.02 | 2.27 | 1.10 | 0.93 | 0.72 |
| LysGly | 0.25 | 2.77 | 21.05 | 1.46 | 1.23 | 0.012 |
| Lys + Ala | 2.14 | 0.85 | 0.74 | 0.46 | 0.38 | 2.89 |
| LysAla | 0.03 | 0.67 | 6.23 | 2.40 | 2.03 | 0.005 |
| Lys + Leu | 1.20 | 1.40 | 1.73 | 0.94 | 0.80 | 0.69 |
| LysLeu | 0.09 | 1.24 | 13.78 | 1.24 | 1.05 | 0.007 |
| Lys + Ser | 0.57 | 0.29 | 0.36 | 0.41 | 0.35 | 1.58 |
| LysSer | 0.19 | 0.48 | 4.16 | 0.29 | 0.25 | 0.046 |
| Lys + Glu | 2.45 | 3.03 | 5.85 | 0.64 | 0.54 | 0.42 |
| LysGlu | 0.24 | 7.77 | 71.46 | 3.52 | 2.97 | 0.003 |
| Lys + Cys | 0.17 | 0.42 | nd | — | — | — |
| LysCys | 0.24 | 1.78 | 5.90 | 2.24 | 1.89 | 0.041 |
| Lys + Phe | 2.17 | 0.90 | 1.16 | 0.32 | 0.27 | 1.87 |
| LysPhe | 0.09 | 1.50 | 11.92 | 2.10 | 1.77 | 0.008 |

Other Table 1 rows (2,3-DMP, TMP, ethyl-, ethenyl-, amino-acid-specific pyrazines) as printed; totals 1.35-111.59;
pyrazines 3.4-87.3 % of total GC-MS area.

### 3.2 Table 2 (MGO pots), selected [M]

2,5(6)-DMP: Lys 0.2 M 3.61, 0.4 M 8.98, LysLys 44.22, Lys+Gly 5.13, LysGly 25.70, Lys+Ala 7.50, LysAla 17.91,
Lys+Leu 6.18, LysLeu 33.44, Lys+Ser 1.71, LysSer 6.22, Lys+Glu 3.96, LysGlu 24.17, Lys+Cys 2.18, LysCys 45.62,
Lys+Phe 6.56, LysPhe 68.22. MPZ 0.03-0.36 (nd in Lys+Glu, Lys+Cys, LysCys); PZ only in Lys+Ser (0.05). TMP
0.15-71.67. 3-ethyl-2,5-DMP 6.77 in Lys + Ala (Strecker acetaldehyde addition) vs 0.62 in LysAla.

### 3.3 Table 3 (glyoxal pots), from the render [M]

| pot | PZ | MPZ | ethylpyrazine |
|---|---:|---:|---:|
| Lys 0.2 M | 6.48 | 0.15 | − |
| Lys 0.4 M | 5.03 | 0.15 | − |
| LysLys | 1.39 | 0.04 | − |
| Lys + Gly | 5.18 | 0.20 | − |
| LysGly | 2.71 | 0.10 | − |
| Lys + Ala | 5.87 | 0.15 | 0.55 |
| LysAla | 1.31 | 0.03 | − |
| Lys + Leu | 1.31 | 0.05 | − |
| LysLeu | 1.19 | 0.06 | − |
| Lys + Ser | 1.82 | 0.05 | − |
| LysSer | 2.17 | 0.05 | − |
| Lys + Glu | 3.02 | 0.06 | − |
| LysGlu | 2.25 | 0.09 | − |
| Lys + Cys | 0.01 | tr | − |
| LysCys | 0.11 | 0.04 | − |
| Lys + Phe | 2.16 | 0.03 | − |
| LysPhe | 2.16 | 0.13 | − |

Recoveries: PZ 0.3, MPZ 0.8, ethylpyrazine 2.3 %. Leu / Phe pots also give 2-(3-methylbutyl)pyrazine (2.10) and
(2-phenylethyl)pyrazine (4.43) with free amino acids only.

### 3.4 Table 4 (Strecker aldehydes, areas x 1e8) [M]

3-methylbutanal: glucose Lys-Leu 9.82 / Lys + Leu 111.98; MGO 0.62 / 81.98; GO 0.02 / 17.75. Phenylacetaldehyde:
glucose Lys-Phe 0.34 / Lys + Phe 4.30; MGO 0.05 / 4.85; GO 0.01 / 0.75.

## 4. Mechanism evidence and its kind

- **Within-study factor** (kind 6), one lab (Ghent): amine form (dipeptide vs free), dicarbonyl identity. Peak areas,
  no replicates of the pots.
- **Fed dicarbonyl** (kind 2): MGO and GO pots give the expected DMP / PZ dominance; amino-acid-specific pyrazines
  need the free amino acid (Strecker aldehyde) [M].
- Scheme 1 (peptide aminoketone without decarboxylation) is [P]; no intermediate detected.
- Lysine alpha vs epsilon reactivity in the dipeptide is hypothesised from Hwang 1994 [C] / [P].

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| EX-052 S1 mixing bound | LEAVES UNCHANGED | all R < 4; non-discriminating by the prereg's rule |
| P2b (extra MPZ route), A-independent flavour | **weakly CONTRADICTS** (~0.5) | R nearly invariant while A/M drops ~10x (sec. 1 answer 2) |
| P1x / P1y vs P2b | **does not distinguish** | |
| CF-077 | LEAVES UNCHANGED | |
| "aminoketone forms only via Strecker" (B18 Strecker step needs a free alpha-amino acid) | CONTRADICTS (direction) | dipeptides with no free alpha-carboxyl make MORE DMP / TMP than free amino acids with MGO and glucose; the engine's Strecker step cannot represent peptide-derived aminoketones |
| evidence map 1.1 (peptides) / 2.7 vanlancker2012 row | SUPPORTS, adds a lab-internal replicate | sec. 6.1 |
| P2a (HCHO / C1 addition to the A-A DHP) | CONTRADICTS a large share in GO + Lys pots | MPZ / PZ ~1-2 % molar |

## 6. Disagreements found

1. **Same nominal pot in Van Lancker 2012 (`vanlancker2012_extraction.md` Table 1)** — Gly + Lys, glucose: PZ / MPZ /
   DMP 1.63 / 2.02 / 2.27 here vs 1.46 / 1.30 / 1.53 (2012); Ala + Lys: 2.14 / 0.85 / 0.74 vs 1.70 / 1.00 / 0.92. R:
   1.10 vs 0.76; 0.46 vs 0.64. Between-paper scatter ~0.15 dex on each area and on R. Same lab, method and
   conditions; cause: run-to-run variation (no replicates printed in either). This is the size of noise to put on
   any Ghent R row.
2. **Van Lancker 2012's dossier lists the charges as unknown** ("only by reference to Van Lancker 2010"). This paper
   gives them: 0.2 M amine / glucose, 0.02 M MGO or GO, 5 mL, initial pH 8, unbuffered. Not a conflict; a gap filled.
3. **R from Ghent (≤ 3.5) vs Leahy lysine-glucose (R ~ 12-17 at 95 °C) and Shibamoto (R ~ 440)**. Candidate causes:
   2,5 / 2,6-DMP co-elution here (inflates DMP, lowers R); pH falling to 4-5 in the unbuffered glucose pots; 130 °C
   vs 75-95 °C; peak areas with theoretical (not measured) recoveries; Leahy's NPD vs EI-MS.
4. **Hwang 1994 as cited** [C]: PZ ~40 % and MPZ ~40 % epsilon-N, 2,6-DMP ≤ 2 %. If 2,5-DMP shared the 2,6 isomer's
   label, the pairing identity label(MPZ) = mean(label(PZ), label(DMP)) would predict ~21 %, not ~40 %. This is a
   second-hand number with 2,5 / 2,6 co-elution; it is a lead for the EX-052 label test (Hwang 1994, JAFC 42:1000, is
   not on disk as a dossier; hwang1995 is).

## 7. Ledger-ready measurements

1. Table 1 PZ / MPZ / 2,5(6)-DMP for the 16 pots with all three detected (peak_area_only, with printed recoveries)
   — mixing-ratio rows with a 100-fold PZ / DMP lever.
2. Table 3 PZ and MPZ in glyoxal-only pots — the cross-product leak (bounds C1-addition routes to MPZ).
3. Table 2 2,5(6)-DMP and TMP, dipeptide vs free (MGO) — the peptide aminoketone channel the engine lacks.
4. Methods charges (0.2 M / 0.02 M, pH 8, 130 °C, 2 h) — to complete the vanlancker2012 rows.
