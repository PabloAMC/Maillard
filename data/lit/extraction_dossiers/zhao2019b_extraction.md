# Zhao, Wang, Xie, Xiao, Cheng, Chen, Wang & Sun 2019 — EXTRACTION
### *"Formation mechanism of aroma compounds in a glutathione-glucose reaction with fat or oxidized fat"*
### Glutathione + glucose 60 + 60 mmol/L in sodium phosphate pH 6.5, ± 53 mg chicken fat or oxidised chicken fat (~1 % w/w), 140 °C, 5 h, sealed 15 mL vial: 49 Maillard / Maillard × lipid volatiles as relative peak areas (MFT and FFT both present), AEDA FD factors, and a 1:1 [13C6]/[12C6]-glucose CAMOLA on 26 compounds in all three systems. One time point, no concentrations, no GSH measurement. Read 2026-10-27.

**Source on disk:** `data/articles/zhao2019b.pdf` (publisher's PDF, 9 pp.). Text from `pdftotext -layout`;
Table 1 checked row by row against a 200 dpi render of p. 438 (matches the text layer); Table 2 from the
text layer (clean). Supplementary material (Table S1, the lipid-only volatiles; Fig. S1, pathways for
Nos. 8, 11, 12, 14, 15) is **not on disk** and was not read.

**Provenance:** **[M]** printed number · **[Q]** qualitative statement · **[F]** read off a figure by me
· **[Z]** derived by me, never printed · **[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | Jian Zhao¹, Tianze Wang¹, Jianchun Xie\*, Qunfei Xiao, Jie Cheng, Feng Chen, Shi Wang, Baoguo Sun (¹ equal contribution) |
| Title | "Formation mechanism of aroma compounds in a glutathione-glucose reaction with fat or oxidized fat" |
| Venue | *Food Chemistry* **270** (2019) **436–444** |
| DOI | 10.1016/j.foodchem.2018.07.106 (printed) |
| Affiliation | Beijing Technology and Business University (BTBU); CAAS Institute of Quality Standard and Testing Technology; Clemson University |
| Dates | received 11 Mar 2018; revised 15 Jul 2018; accepted 17 Jul 2018; online 18 Jul 2018 |

Same group as Wang T. et al. 2020 (`wang2020_extraction.md`, which reuses these conditions for the
intermediates) and close to Du et al. 2025 (`du2025_extraction.md`, CAAS, GSH + glucose + free C18 fatty
acids). Cited in `sun2022_extraction.md` pointer table as item 2.

---

## 1. THE ONE-PARAGRAPH ANSWER

0.3 mmol GSH + 0.3 mmol glucose in 5 mL sodium phosphate pH 6.5 (**60 + 60 mmol/L** [Z]), heated 5 h at
140 °C in a sealed vial, with or without 53 mg chicken fat or oxidised chicken fat (PV 237 meq O₂/kg).
Results are **relative peak areas** (compound / an internal standard injected into the GC, not into the
sample), so only ratios between systems are usable. **Fat lowers both meaty thiols: MFT 0.52 → 0.32
(fat, −38 %) → 0.18 (oxidised fat, −65 %); FFT 0.88 → 0.63 (−28 %) → 0.27 (−69 %)** [Z], all letters
distinct. In the same direction, most glucose-derived sulfur heterocycles fall (thieno[3,2-b]thiophene
5.19 → 1.17 → 0.88), while **lipid-only-carbon alkylthiophenes appear** (2-pentylthiophene 0 → 1.04 →
8.91; 2-butyl-, 2-hexyl-, 2-heptyl-, 2-propylthiophene; 1-hexanethiol only with oxidised fat), i.e. lipid
carbonyls divert H₂S. CAMOLA: **MFT and FFT are ~100 % a glucose C5 fragment in all three systems**
(unlabelled : 13C5 ≈ 50 : 47); fat changes the isotope pattern of none of the Maillard products except
2-ethylthiophene (glucose share ~90 % → ~44 % with oxidised fat). Thiazoles carry **non-glucose C3**
(2-acetylthiazole 74 % unlabelled vs 18 % 13C3), which the authors attribute to MGO/hydroxyacetone from
GSH's glutamyl residue. **This sign (fat lowers MFT) is the opposite of du2025** (free C18 acids raised MFT
by 19–30 % at 145 °C, 2 h).

---

## 2. METHODS [M] (pp. 437–438)

| item | as printed |
|---|---|
| Reactants | glutathione (99 %), [12C6]-D-glucose (98 %), [13C6]-D-glucose (99 %), J&K; refined odourless chicken fat (Tianjin Mu Yang Oil Factory) |
| Oxidised fat | prepared from the same chicken fat "according to Yang et al. (2015)"; **peroxide value 237 meq O₂/kg, acid value 2.00 mg KOH/g**; oxidation conditions not restated [!] |
| Systems | (a) GSH-Glucose (control); (b) GSH-Glucose-Chicken fat; (c) GSH-Glucose-Oxidised chicken fat; (d) chicken fat alone; (e) oxidised chicken fat alone (blanks; their volatiles are in Table S1, not on disk) |
| Charge | **GSH 0.3 mmol + glucose 0.3 mmol in 5 mL sodium phosphate buffer (pH 6.5)**; fat **53.0 mg** ("ca. 1 % w/w of the total content") [Z: 60 mmol/L each; GSH 18.4 g/L] |
| Buffer molarity | **not stated** [!]. 60 mmol/L GSH (an acid) in a buffer of unknown strength very likely lowers the start pH below 6.5; no measured pH, before or after [!] |
| Vessel | **15 mL Synthware glass vial, screw cap** (headspace ≈ 10 mL air; not flushed) |
| Heating | Parallel Synthesis Poly-block 4 (HEL), **140 °C, with stirring, 5 h**; come-up not stated; cooled in cold water before opening |
| Replicates | triplicate reactions; "results were the means of three replicates" |
| SPME | 50/30 µm CAR/PDMS/DVB; sample equilibrated **50 °C, 10 min**, fibre in headspace **20 min**, stirred |
| GC-MS | 7890B/5975C, DB-Wax 30 m; 40 → 70 °C at 5 °C/min, → 160 °C at 2 °C/min, → 230 °C at 8 °C/min |
| "Internal standard" | **0.1 µL of 1,2-dichlorobenzene (60 ng/µL in methanol) injected into the GC inlet immediately before fibre desorption** — i.e. an injection standard, not a sample internal standard. It corrects nothing about extraction [!] |
| Quantity reported | **relative peak area = compound area / IS area** (dimensionless) |
| GC-O / AEDA | DATU 2000 olfactometer, three sniffers; AEDA on system (a) only, split 1:2 … 1:128 then shortened fibre (FD 1–512) |
| CAMOLA | 0.15 mmol [13C6]-glucose + 0.15 mmol [12C6]-glucose replacing 0.3 mmol glucose, in (a), (b), (c); isotopomer areas M⁺ … M⁺+n normalised, corrected for 1.10 % 13C, 0.76 % 33S, 4.20 % 34S and for (M⁺−1)/M⁺ H-loss |
| Statistics | one-way ANOVA, Duncan, p < 0.05 (SPSS 19) |

---

## 3. RESULTS

### 3.1 Table 1 (p. 438): volatiles, relative peak area (mean ± SD, n = 3), letters as printed [M]

a = GSH-Glucose; b = + chicken fat; c = + oxidised chicken fat. "–" = not detected. FD from AEDA on (a).
ID: RI, MS (NIST 15 + interpretation), O (odour), S (authentic standard).

| RI | Compound | a | b | c | Odour | FD | ID |
|---|---|---|---|---|---|---|---|
| 1024 | thiophene | – | – | – | garlic | 256 | RI, O |
| 1109 | 2-methylthiophene | 6.77 ± 0.56a | 5.00 ± 0.31b | 2.70 ± 0.12c | sulfury | 256 | RI, MS, O, S |
| 1160 | 2,5-dimethylthiophene | 1.29 ± 0.16a | 0.60 ± 0.08b | 0.26 ± 0.02c | sulfury | 1 | RI, MS, O, S |
| 1174 | 2-ethylthiophene | 0.59 ± 0.11a | 0.35 ± 0.01b | 0.52 ± 0.02a | meaty | 32 | RI, MS, O, S |
| 1149 | 1-hexanethiol | – | – | 0.64 ± 0.05 | green, meaty | – | RI, MS, O |
| 1225 | thiazole | – | – | – | meaty | 64 | RI, O |
| 1214 | 2,3-dimethylthiophene | 0.31 ± 0.03a | 0.17 ± 0.02b | – | meaty | 128 | RI, MS, O |
| 1249 | 2-propylthiophene | – | 0.63 ± 0.02b | 1.28 ± 0.01a | green, meaty | – | RI, MS, O, S |
| 1302 | **2-methyl-3-furanthiol** | **0.52 ± 0.07a** | **0.32 ± 0.02b** | **0.18 ± 0.02c** | meaty | 128 | RI, MS, O, S |
| 1311 | 2,4-dimethylthiazole | 0.17 ± 0.01a | 0.14 ± 0.01b | – | meaty | 16 | RI, MS, O, S |
| 1333 | 3-mercapto-2-pentanone | – | – | – | meaty | 512 | RI, O, S |
| 1347 | 2-butylthiophene | – | 0.56 ± 0.02b | 2.50 ± 0.27a | green, meaty | – | RI, MS, O, S |
| 1356 | 2-ethyl-4-methylthiazole | 0.14 ± 0.02 | – | – | – | – | RI, MS |
| 1371 | 2,4,5-trimethylthiazole | 0.60 ± 0.04a | 0.30 ± 0.04b | 0.16 ± 0.03c | earthy, sulfury | 4 | RI, MS, O, S |
| 1419 | **2-furfurylthiol** | **0.88 ± 0.05a** | **0.63 ± 0.07b** | **0.27 ± 0.03c** | coffee | 512 | RI, MS, O, S |
| 1445 | furfural | 0.41 ± 0.03a | 0.46 ± 0.01a | 0.45 ± 0.02a | caramel | 4 | RI, MS, O, S |
| 1446 | 2-pentylthiophene | – | 1.04 ± 0.17b | 8.91 ± 0.78a | green, meaty | – | RI, MS, O, S |
| 1472 | tetrahydrothieno[3,2-b]furan | 0.25 ± 0.06a | 0.09 ± 0.01b | – | rice | 256 | MS |
| 1485 | 2-acetylfuran | 0.33 ± 0.01c | 0.57 ± 0.03b | 0.87 ± 0.06a | – | – | RI, MS |
| 1505 | 2-methyltetrahydrothiophen-3-one | 0.03 ± 0.00c | 0.06 ± 0.01b | 0.09 ± 0.01a | meaty | 4 | RI, MS, O, S |
| 1535 | tetrahydrothiophen-3-one | 0.18 ± 0.04a | 0.15 ± 0.02a | 0.14 ± 0.02a | meaty | 8 | RI, MS, O, S |
| 1538 | 2-methyl-5-[(methylthio)methyl]furan | 0.25 ± 0.03a | 0.09 ± 0.01b | – | meaty | 64 | RI, MS, O |
| 1549 | 2-hexylthiophene | – | 1.27 ± 0.09b | 3.33 ± 0.25a | green, meaty | – | RI, MS, O, S |
| 1565 | 2-pentylpyridine | – | 0.30 ± 0.05b | 3.99 ± 0.09a | fatty, roasted | – | RI, MS, O, S |
| 1573 | anti-3,5-dimethyl-1,2,4-trithiolane | 2.91 ± 0.38a | 1.95 ± 0.21b | 1.58 ± 0.08b | garlic | 32 | RI, MS, O |
| 1552 | 5-methylfurfural | 3.33 ± 0.17c | 3.93 ± 0.22b | 4.90 ± 0.14a | caramel | 128 | RI, MS, O, S |
| 1594 | syn-3,5-dimethyl-1,2,4-trithiolane | 2.75 ± 0.21a | 1.73 ± 0.17b | 1.47 ± 0.04b | garlic | 32 | RI, MS, O |
| 1622 | 2-acetylthiazole | 0.69 ± 0.02a | 0.73 ± 0.01a | 0.63 ± 0.04b | roasted, nutty | 32 | RI, MS, O, S |
| 1647 | furfuryl alcohol | 0.32 ± 0.03b | 0.44 ± 0.07b | 0.60 ± 0.04a | caramel | 16 | RI, MS, O, S |
| 1655 | 2-heptylthiophene | – | 0.48 ± 0.05b | 3.34 ± 0.13a | fatty, meaty | – | RI, MS, O |
| 1664 | 2-(2-furanylmethyl)-5-methylfuran | 0.53 ± 0.01 | – | – | – | – | RI, MS |
| 1662 | thiophene-2-carboxaldehyde | 1.39 ± 0.08a | 0.75 ± 0.05b | 0.54 ± 0.07c | meaty | 16 | RI, MS, O, S |
| 1738 | 3-acetylthiophene | 0.59 ± 0.02a | 0.59 ± 0.04a | 0.57 ± 0.02a | meaty | 16 | RI, MS, O, S |
| 1744 | 2-acetylthiophene | 0.64 ± 0.04a | 0.63 ± 0.03a | 0.59 ± 0.01a | meaty | 16 | RI, MS, O, S |
| 1745 | 5-methylthiophene-2-carboxaldehyde | 0.98 ± 0.04a | 0.79 ± 0.06b | 0.30 ± 0.01c | meaty | 32 | RI, MS, O, S |
| 1777 | 3-methylthiophene-2-carboxaldehyde | 1.55 ± 0.06a | 1.42 ± 0.06b | 1.33 ± 0.12b | meaty | 32 | RI, MS, O, S |
| 1790 | 2-acetyl-3-methylthiophene | 0.27 ± 0.05a | 0.22 ± 0.02b | 0.09 ± 0.00c | meaty | 64 | RI, MS, O |
| 1817 | 1-(2-thienyl)-1-propanone | – | – | 0.16 ± 0.02 | burnt, fatty | – | RI, MS, O |
| 1837 | thieno[3,2-b]thiophene | 5.19 ± 0.21a | 1.17 ± 0.11b | 0.88 ± 0.09b | meaty | 64 | RI, MS, O |
| 1856 | 2-mercaptophenol | 0.56 ± 0.10a | 0.26 ± 0.05b | – | meaty | 32 | MS, O |
| 1868 | 1(3H)-isobenzofuranthione | 1.21 ± 0.14a | 0.36 ± 0.02b | 0.33 ± 0.02b | meaty | 32 | MS, O |
| 1901 | 2-methylthieno[3,2-b]thiophene | 4.91 ± 0.40a | 0.39 ± 0.01b | 0.34 ± 0.02b | meaty | 32 | RI, MS, O |
| 1918 | 2-thiophenemethanol | 0.50 ± 0.06a | 0.51 ± 0.04a | 0.39 ± 0.04b | meaty, burnt | 8 | RI, MS, O |
| 1921 | 2-methylthieno[2,3-b]thiophene | 0.80 ± 0.14a | 0.13 ± 0.01b | 0.09 ± 0.01b | meaty | 256 | RI, MS, O |
| 1932 | 2H-thiopyran[2,3-b]-2H-thiete | 0.87 ± 0.16a | 0.45 ± 0.01b | 0.28 ± 0.03c | meaty | 512 | MS, O |
| 1943 | 2-acetylpyrrole | 0.91 ± 0.14a | 0.48 ± 0.01b | 0.51 ± 0.05b | – | – | RI, MS |
| 1986 | 1-(2-thienyl)-1-pentanone | – | – | 0.57 ± 0.07 | – | – | RI, MS |
| 1991 | 2-ethylthieno[3,2-b]thiophene | 1.77 ± 0.34a | 0.06 ± 0.02b | 0.04 ± 0.00b | meaty | 8 | RI, MS, O |
| 2072 | 2,2′-bithiophene | 0.58 ± 0.10 | – | – | meaty | 64 | MS, O |
| 2094 | 2-formyl-2,3-dihydrothiophene | 0.94 ± 0.08a | 0.76 ± 0.04b | 0.72 ± 0.08b | burnt, meaty | 512 | MS, O |
| 2096 | 1-(2-thienyl)-1-hexanone | – | – | 0.36 ± 0.04 | – | – | RI, MS |
| 2140 | 2,3′-bithiophene | 0.72 ± 0.04a | 0.12 ± 0.00b | – | meaty | 32 | MS, O |
| 2249 | 2,3-dihydro-3,5-dihydroxy-6-methyl-4H-pyran-4-one | 0.15 ± 0.02a | 0.14 ± 0.05a | 0.07 ± 0.05b | – | – | RI, MS |

Notes [!]: (i) thiophene, thiazole and 3-mercapto-2-pentanone (FD 256, 64, **512**) were smelled but
never seen by MS ("–" in all three columns), so they have no amount. (ii) The RI column is not monotone
(1149 after 1174; 1552 after 1573; 1664 after 1655; 1662 after 1664) — row order is the authors', not
elution order. (iii) Text: "thirty-three sulfur-containing compounds" in (a); the table supports that
count only if GC-O-only rows are included [Z].

**Fat effect on the meaty thiols and key S-heterocycles [Z]** (ratio to system a):

| compound | b/a (chicken fat) | c/a (oxidised fat) |
|---|---|---|
| 2-methyl-3-furanthiol | 0.62 | 0.35 |
| 2-furfurylthiol | 0.72 | 0.31 |
| 2-methylthiophene | 0.74 | 0.40 |
| thiophene-2-carboxaldehyde | 0.54 | 0.39 |
| thieno[3,2-b]thiophene | 0.23 | 0.17 |
| 2-methylthieno[3,2-b]thiophene | 0.08 | 0.07 |
| anti + syn-3,5-dimethyl-1,2,4-trithiolane | 0.65 | 0.54 |
| 2-acetylthiazole | 1.06 | 0.91 |
| furfural / 5-methylfurfural / 2-acetylfuran (non-S) | 1.12 / 1.18 / 1.73 | 1.10 / 1.47 / 2.64 |

Σ lipid-derived 2-alkylthiophenes + 1-hexanethiol + acylthiophenes (rows absent from a) [Z]: b 3.98, c
20.59 relative area units; compare MFT + FFT: a 1.40, b 0.95, c 0.45.

### 3.2 Table 2 (p. 441): CAMOLA isotopomer distribution, % of the molecular-ion cluster [M]

Glucose 1:1 [13C6]/[12C6]. Column = number of 13C atoms. Shown for the compounds with 13C isotopomers and
FD ≥ 16; data for systems b and c printed **only for 2-ethylthiophene** (the text says all other patterns
were unchanged by fat; numbers "not shown") [!].

| No. | Compound | M⁺ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2-acetylthiophene | 126 | 51 | <1 | <1 | 2 | 3 | <1 | 43 | | |
| 2 | 5-methylfurfural | 110 | 49 | <1 | <1 | <1 | 2 | <1 | 47 | | |
| 3 | 5-methylthiophene-2-carboxaldehyde | 126 | 52 | <1 | <1 | 1 | 2 | 3 | 42 | | |
| 4 | **2-methyl-3-furanthiol** | 114 | **50** | <1 | <1 | <1 | 2 | **47** | | | |
| 5 | **2-furfurylthiol** | 114 | **52** | <1 | <1 | <1 | <1 | **47** | | | |
| 6 | thiophene-2-carboxaldehyde | 112 | 51 | 2 | 3 | 2 | 2 | 40 | | | |
| 7 | furfuryl alcohol | 98 | 50 | 4 | 4 | 1 | <1 | 41 | | | |
| 8 | 2-formyl-2,3-dihydrothiophene | 114 | 50 | <1 | <1 | 3 | 3 | 43 | | | |
| 9 | 2-methylthiophene | 98 | 48 | 1 | 3 | 3 | 2 | 41 | | | |
| 10 | 2-methyl-5-[(methylthio)methyl]furan | 142 | 51 | <1 | <1 | <1 | <1 | 2 | 46 | <1 | |
| 11 | 2-acetyl-3-methylthiophene | 140 | 49 | 3 | <1 | <1 | <1 | 44 | <1 | | |
| 12 | 3-acetylthiophene | 126 | 49 | 1 | 2 | 3 | 42 | 1 | 2 | | |
| 13 | 3-methylthiophene-2-carboxaldehyde | 126 | 51 | <1 | 1 | 3 | 42 | 2 | <1 | | |
| 14 | 2,3-dimethylthiophene | 112 | 48 | 4 | 4 | 1 | 43 | <1 | <1 | | |
| 15a | 2-ethylthiophene (system a) | 112 | 52 | 1 | 1 | <1 | 5 | 1 | 40 | | |
| 15b | 2-ethylthiophene (system b) | 112 | 59 | 1 | <1 | <1 | 5 | <1 | 33 | | |
| 15c | 2-ethylthiophene (system c) | 112 | 75 | 1 | <1 | <1 | 3 | <1 | 19 | | |
| 16 | tetrahydrothieno[3,2-b]furan | 128 | 50 | 2 | <1 | <1 | <1 | 5 | 42 | | |
| 17 | 2,2′-bithiophene | 166 | 46 | 2 | 4 | 1 | 1 | 3 | 43 | <1 | |
| 18 | 2,3′-bithiophene | 166 | 48 | <1 | 4 | <1 | 2 | 2 | 44 | <1 | |
| 19 | 2-methylthieno[3,2-b]thiophene | 154 | 48 | 1 | 2 | <1 | 1 | 44 | 2 | 1 | |
| 20 | 2-methylthieno[2,3-b]thiophene | 154 | 47 | 5 | <1 | <1 | <1 | 44 | <1 | 3 | |
| 21 | thieno[3,2-b]thiophene | 140 | 54 | <1 | 3 | <1 | 9 | 4 | 28 | | |
| 22 | 2H-thiopyran[2,3-b]-2H-thiete | 140 | 50 | <1 | 4 | 3 | 27 | 4 | 10 | | |
| 23 | 2,4-dimethylthiazole | 113 | 66 | 2 | 3 | 24 | 1 | 4 | | | |
| 24 | 2-acetylthiazole | 127 | 74 | 1 | 1 | 18 | 3 | 2 | | | |
| 25 | 2-mercaptophenol | 126 | 48 | <1 | 3 | 2 | 3 | 3 | 41 | | |
| 26 | 1(3H)-isobenzofuranthione | 150 | 51 | <1 | <1 | <1 | <1 | 3 | 45 | | |

Text [Q]: the ten lipid-borne alkyl compounds of Table 1 (1-hexanethiol, 2-propyl-/butyl-/pentyl-/hexyl-/
heptylthiophene, three 1-(2-thienyl)-1-alkanones, 2-pentylpyridine) and **both 3,5-dimethyl-1,2,4-
trithiolanes showed no 13C isotopomers** ("data not shown").

**Derived carbon origin [Z]** (glucose share ≈ 2 × labelled isotopomers, 1:1 mix, equal reactivity assumed):

| compound | glucose-carbon share | non-glucose share | reading |
|---|---|---|---|
| MFT (C5) | ~94–98 % (13C5, + 2 % 13C4) | ≤ 6 % | glucose C5 fragment (authors: 1-deoxypentosone from glucose) |
| FFT (C5) | ~94 % (13C5) | ~6 % | glucose C5 fragment (authors: 3-deoxypentosone from glucose) |
| 2-ethylthiophene, a / b / c | ~90 / ~76 / ~44 % (13C6 + 13C4) | ~10 / ~24 / ~56 % | lipid carbon enters only with fat (authors: 2,4-hexadienal + H₂S) |
| 2-acetylthiazole (C3 module) | ~36 % | **~64 %** | non-glucose C3 (authors: MGO/hydroxyacetone from GSH's Glu) |
| 2,4-dimethylthiazole | ~48 % | ~52 % | same |
| thieno[3,2-b]thiophene | ~74 % (13C6 + 13C4) | ~26 % | C4 + C2 route present |
| 3,5-dimethyl-1,2,4-trithiolanes | 0 % | 100 % | acetaldehyde not from glucose |

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Model context.** `GSH` is a sulfide source only: `r_gsh_h2s`, GSH → H₂S + FRAG_C(10) + FRAG_N(3),
`k_gsh_h2s` from Zheng & Ho 1994 Table III, pH-resolved (`src/kinetic_core/sulfur.py` l. 246,
`parameters_sulfur.py` l. 1934). Glucose is the hexose route. The lipid lane makes oxidation aldehydes;
no aldehyde + H₂S → alkylthiophene sink exists. Nothing here was run against the model.

**(1) Direction test, lipid × sulfur (usable).** In a GSH + glucose pot, **triglyceride fat lowers MFT
(−38 %) and FFT (−28 %); oxidised fat lowers them further (−65 %, −69 %)**. The ratios are robust to the
IS problem only insofar as matrix effects on SPME are similar across systems (53 mg fat in 5 mL changes
the headspace partition of lipophilic thiols — itself a candidate cause of the drop [!]).
**Conflict:** du2025 (same chemistry family, free C18 fatty acids at 1 : 20 FA : GSH, 145 °C, 2 h) found
MFT **up** 19–30 %. Differences: free acid vs triglyceride; 0.05 mmol FA vs 53 mg fat (~0.06 mmol
triglyceride, ~0.18 mmol acyl chains); 2 h vs 5 h; oxidised fat (PV 237) here. Read together: the sign of
the lipid effect on MFT is **not established**; the oxidised-fat (pre-formed carbonyl) arm is the one
that consistently suppresses sulfur volatiles. The same suppression appears in Cys + xylose with
phospholipids (`zhang2021_extraction.md`: MFT −40 to −64 %).

**(2) Sink evidence (qualitative to semi-quantitative).** Alkylthiophenes with no glucose carbon rise to
~20 relative-area units with oxidised fat, ~15 × the MFT + FFT left in that pot. If an aldehyde + H₂S →
2-alkylthiophene sink is ever added, this pot gives a *direction and order-of-magnitude ratio* check
(sink products ≫ meaty thiols under oxidised fat), not a rate.

**(3) Carbon origin (usable, qualitative).** MFT and FFT from GSH + glucose are glucose C5 fragments
(~95 %); fat does not open a second MFT/FFT carbon route. This supports routing hexose → C5 deoxyosone →
MFT/FFT in the model and contradicts nothing in the current network. Note the contrast with
`wang2012_extraction.md` (GSH + xylose + thiamine: MFT > 99 % non-xylose, i.e. thiamine-borne) — the GSH
lane does not supply MFT carbon; the sugar or thiamine does.

**(4) Glutathione's non-H₂S chemistry (qualitative).** Two observations argue that GSH contributes more
than H₂S: (i) thiazoles carry ~50–64 % non-glucose C3, assigned (speculatively) to Glu-derived MGO;
(ii) the dimethyltrithiolanes are 100 % non-glucose and present **without fat**, so their acetaldehyde
cannot be lipid-derived as the authors write [!] — GSH (Cys Strecker/thermolysis) is the only other
carbon source. Neither is a rate.

**(5) Not usable.** Any absolute amount (no sample internal standard; relative peak areas only); any
yield; anything about GSH consumption or Amadori formation (not measured); pH (not measured); the
proposed pathways (Figs. 1–4) beyond what Table 2 shows; the lipid-only blanks (Table S1 not on disk).

**Caveats, collected.**
1. Single time point (5 h, 140 °C) — far past the GSH half-life implied by du2025/wang2020/zhou2023c.
2. Buffer molarity and actual pH unknown; 60 mmol/L GSH acidifies.
3. "Internal standard" injected into the inlet, not the sample: values are area ratios, not
   concentrations; fat changes the SPME matrix.
4. Oxidised-fat preparation only by reference (Yang et al. 2015); PV/AV printed.
5. CAMOLA numbers for fat systems printed only for 2-ethylthiophene.
6. Three GC-O-only odorants (incl. FD-512 3-mercapto-2-pentanone) have no MS amount.
