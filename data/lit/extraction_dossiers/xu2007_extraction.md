# Xu, Gao, Liu & Zhao 2008 (online 2007) — EXTRACTION
### *"Effects of supercritical carbon dioxide on volatile formation from Maillard reaction between ribose and cysteine"*
### Ribose 50 mM + cysteine 50 mM in 0.2 M "phosphate" buffer (made from sodium **pyrophosphate** salts) at initial pH **5.6, 7.0, 8.0**, 300 mL stainless-steel vessel, **140 °C / 1 h**, either conventional ("Contr.") or under SC-CO2 at 40 MPa; HS-SPME GC-MS, ~70 volatiles in **ng/mL headspace-equivalent of mixture** (tridecane, RF = 1). **MFT 16.62 / 2.63 / – and FFT 22.09 / – / – ng/mL (control, initial pH 5.6 / 7.0 / 8.0; final 5.5 / 6.2 / 6.8).** The pH 5.6 columns are the same data as `xu2008`'s 1:1 columns. Read 2026-10-28.

**Full read of `data/articles/xu2007.pdf`** (8 pp., born-digital Wiley PDF). Prose and Table 1 from
`pdftotext -layout`; the Table 1 thiol and thienothiophene rows were spot-checked against a 200 dpi render of
p. 332 and agree. Fig. 1 (apparatus) and Fig. 2 (A280 / A420 vs pH, error bars, n = 2) were not digitised
[!].

**Year:** printed "J Sci Food Agric **88**:328–335 (**2008**)", "Published online 1 October 2007", © 2007 SCI.
Crossref: issue 2, 2008-01-30. The repo stem `xu2007` reflects the online date; cite as **Xu et al. 2008**
(JSFA 88:328). `sun2022_extraction.md` records a citing review that misprints it as "2010".

**Provenance:** **[M]** printed · **[Q]** qualitative · **[Z]** derived · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value |
|---|---|
| Authors | Honggao Xu, Yanxiang Gao* (China Agricultural University, Beijing), Xuan Liu (same), Jian Zhao (Charles Sturt University, Wagga Wagga) |
| Title | "Effects of supercritical carbon dioxide on volatile formation from Maillard reaction between ribose and cysteine" |
| Venue | *J. Sci. Food Agric.* **88**: 328–335 (2008); issue 2 |
| DOI | **10.1002/jsfa.3093** (printed) |
| Dates | received 24 Feb 2007; revised 2 Jul 2007; accepted 23 Jul 2007; online 1 Oct 2007 |
| Grant | NSFC 20436020 |

---

## 1. THE ONE-PARAGRAPH ANSWER

The pH-axis sibling of `xu2008` (FRI 41:730, ratio axis) and the parent of `xu2010` (pressure axis). In the
**conventional** 140 °C / 1 h pot, MFT falls **16.62 → 2.63 → below detection** and FFT **22.09 → – → –**
ng/mL as the initial pH goes 5.6 → 7.0 → 8.0 (final 5.5 / 6.2 / 6.8), while pyrazines, 2-acetylthiazole and
3,5-dimethyl-1,2,4-trithiolanes rise. A dihydrothienothiophene (LRI 1325) dominates every headspace (441 /
683 / 567). SC-CO2 (40 MPa) acidifies the pot (final 5.5 / 6.0 / 6.3), suppresses most volatiles but
**raises the thiols** (MFT 34.4 / 25.8 / 7.1; FFT 25.1 / 14.6 / 1.8) and furfural, and eliminates the
trithiolanes. The authors call the quantities approximate headspace amounts, "not … absolute concentrations
in the aqueous solutions."

---

## 2. METHODS [M]

| item | as printed |
|---|---|
| Buffer | "0.2 mol L⁻¹ phosphate buffer" prepared "by mixing appropriate volumes of 0.2 mol L⁻¹ **sodium pyrophosphate dibasic** and 0.2 mol L⁻¹ **sodium pyrophosphate tetrabasic**" — i.e. **pyrophosphate**, as `xu2008` states explicitly [!] |
| Stocks | ribose **7.5064 g** or cysteine **6.0576 g** in ~100 mL buffer, pH re-adjusted with tetrabasic pyrophosphate, made up to **500 mL** → 0.100 M each [Z ✓: 7.5064/150.13/0.5 = 0.1000; 6.0576/121.16/0.5 = 0.1000] |
| Pot | **50 mL ribose stock + 50 mL cysteine stock** mixed in the 300 mL reactor → **50 mM + 50 mM in 100 mL** (Table 1 footnote b: "5 mmol ribose and 5 mmol cysteine in 100 mL 0.2 mol L⁻¹ phosphate buffer") |
| SC-CO2 arm | CWYF-2 apparatus; CO2 through activated carbon; **40 MPa, 140 °C, held 1 h**; cooled with −3 °C 30 % ethanol; depressurised through **three 50 mL buffer absorbers**; volatiles = reaction mixture + absorbers combined |
| Control arm | "heated at 140 °C for the same period using the same reaction apparatus … but without SC-CO2" (pressure/headspace gas not stated [!]) |
| Replicates | treatments in duplicate; analyses ≥ triplicate; "mean CV … < 20 %"; no compound CV > 47 % except small peaks |
| HS-SPME | 1 mL mixture (15 mL absorber) in 15 mL (60 mL) bottle + **15.13 ng tridecane**; 10 min at 60 °C, then DVB/CAR/PDMS 50/30 µm 20 min at 60 °C |
| GC-MS | Agilent 6890N/5973, DB-5ms 60 m × 0.32 mm × 0.5 µm; 40 → 60 °C (20 °C/min, 5 min) → 4 °C/min → 250 °C 10 min; splitless 270 °C 10 min |
| Quantity | TIC area vs tridecane, **response factor 1**; "Approximate quantities in headspace (**ng mL⁻¹ of mixture**)"; tr < 0.3, "−" < ~0.07 ng/mL |
| Browning | A280, A420 (Fig. 2) — not digitised |

---

## 3. TABLE 1 [M] — "Approximate quantities of volatiles identified in the headspace of ribose–cysteine model system under different conditions" (ng mL⁻¹ of mixture)

| LRI | Compound | Acid Contr. | Acid SC-CO2 | Neutral Contr. | Neutral SC-CO2 | Alkaline Contr. | Alkaline SC-CO2 | ID |
|---|---|---|---|---|---|---|---|---|
| | **Initial pH** | 5.6 | 5.6 | 7.0 | 7.0 | 8.0 | 8.0 | |
| | **Final pH** | **5.5** | **5.5** | **6.2** | **6.0** | **6.8** | **6.3** | |
| | *Thiophene* | | | | | | | |
| 773 | 2-Methylthiophene | 12.64 | − | 4.28 | − | 0.54 | − | ms + LRI |
| 873 | 2,5-Dimethylthiophene | 9.74 | − | 6.68 | − | − | − | ms + LRI |
| 891 | 2,3-Dimethylthiophene | 0.95 | − | 1.46 | − | − | − | ms + LRI |
| 1068 | 2-Thienylmethanol | 22.96 | 19.75 | 8.43 | 25.13 | 4.71 | 30.12 | ms |
| 1095 | 2-Acetylthiophene | 3.78 | 0.70 | 9.34 | 1.09 | 11.78 | 1.74 | ms + LRI |
| 1100 | 2-Formyl-5-methyl-2,5-dihydrothiophene | 8.25 | 2.66 | 4.22 | 5.35 | 3.65 | 6.73 | ms |
| 1128 | 5-Methyl-2-thiophenecarboxaldehyde | 28.05 | 19.36 | 12.38 | 38.45 | 4.55 | 34.33 | ms + LRI |
| 1149 | 2-Butylthiophene | 3.13 | − | 3.6 | − | 11.64 | 0.62 | ms |
| 1162 | 2-Acetyl-5-methylthiophene | − | − | − | − | − | 0.88 | ms + LRI |
| 1176 | 1-(2-Thienyl)-1-butanone | + | − | 4.58 | tr | 26.73 | 1.85 | ms |
| 1189 | 2-Propionylthiophene | 15.23 | 5.31 | 7.5 | 12.06 | 20.82 | 16.46 | ms + LRI |
| 1194 | 2-Acetyl-3-methylthiophene | + | − | 12 | tr | 18.64 | 1.91 | ms |
| 1208 | 3-Ethyl-2-formylthiophene | 8.9 | 10.88 | 2.98 | 11.04 | 2.23 | 10.46 | ms + LRI |
| 1249 | Dimethylformylthiophene | 3.66 | 6.22 | 3.27 | 5.1 | 4.64 | 3.79 | ms |
| 1278 | 2-Methoxy-3,4,5-trimethylthiophene | + | 1.04 | 2.84 | 2.91 | 3.61 | 2.77 | ms |
| 1289 | 3-Acetyl-2,5-dimethylthiophene | 6.9 | − | 21.49 | 1.58 | 13.75 | 2.49 | ms |
| | *Thiophenone* | | | | | | | |
| 959 | Dihydro-3-(2H)-thiophenone | 1.58 | − | 1.42 | − | 1.33 | 0.94 | MS + LRI |
| 996 | Dihydro-2-methyl-3(2H)-thiophenone | 31.75 | 11.38 | 11.08 | 13.96 | 7.18 | 19.52 | MS + LRI |
| 1026 | 4,5-Dihydro-2,4-dimethyl-3(2H)-thiophenone (E or Z) | 6.97 | 0.32 | 2.46 | 1.13 | − | 1.71 | ms + LRI |
| | *Fused bicyclic* | | | | | | | |
| 1197 | 2,3-Dihydro-6-methylthieno[2,3-c]furan | 56.08 | 13.13 | 99.29 | 38.09 | 52.45 | 42.93 | ms + LRI |
| 1220 | Thieno[2,3-b]thiophene | 17.36 | 9.04 | 25.85 | 7.2 | 14.1 | 8.26 | ms + LRI |
| 1262 | Thieno[3,2-b]thiophene | 10.26 | 1.36 | 29.26 | 3.78 | 34.78 | 4.25 | ms |
| 1306 | 2-Methyl-thieno[2,3-b]thiophene | 2.46 | − | 2.46 | − | 1.72 | − | ms |
| 1311 | a-Dihydrothienothiophene | 11.39 | − | 9.48 | 2.22 | 6.69 | 4.17 | ms + LRI |
| 1325 | a-Dihydrothienothiophene | 441.52 | 62.25 | 683.00 | 180.70 | 567.31 | 294.30 | ms |
| 1364 | 2-Methylthieno[3,2-b]thiophene | 21.20 | 5.87 | 39.43 | 11.05 | 110.41 | 13.72 | ms |
| 1381 | a-Methyldihydrothienothiophene | 75.56 | 4.26 | 58.52 | 12.12 | 32.96 | 27.1 | ms + LRI |
| 1418 | a-Methyldihydrothienothiophene | 29.57 | 9.95 | 7.26 | 23.29 | 9.52 | 32.69 | ms + LRI |
| 1423 | a-Methyldihydrothienothiophene | 39.04 | 25.27 | 33.15 | 34.8 | 30.73 | 49.67 | ms + LRI |
| 1477 | a-Dimethyldihydrothienothiophene | 36.62 | 5.52 | 19.34 | 11.58 | 12.87 | 17.23 | ms |
| | *Thiol* | | | | | | | |
| 869 | **2-Methyl-3-furanthiol** | **16.62** | **34.42** | **2.63** | **25.84** | **−** | **7.06** | ms + LRI |
| 913 | **2-Furanmethanethiol** | **22.09** | **25.13** | **−** | **14.57** | **−** | **1.79** | MS + LRI |
| 980 | 3-Thiophenethiol | 25.55 | 36.94 | 2.19 | 19.12 | − | 3.63 | ms + LRI |
| 1066 | 2-Methyl-3-thiophenethiol | 25.37 | 14.22 | − | 9.93 | − | 2.69 | ms + LRI |
| | *Polysulfur heterocyclic* | | | | | | | |
| 1153 | 3,5-Dimethyl-1,2,4-trithiolane (E or Z) | 5.02 | − | 18.86 | − | 33.31 | − | MS |
| 1160 | 3,5-Dimethyl-1,2,4-trithiolane (E or Z) | 5.90 | − | 19.7 | − | 31.73 | − | MS |
| 1185 | 1,2-Dithian-4-one | 8.1 | 3.37 | 1.92 | 6.46 | − | 9.75 | ms |
| 1232 | 3-Methyl-1,2-dithian-4-one | 8.97 | 40.12 | 1.47 | 51.19 | 1.31 | 65.59 | MS |
| 1264 | 3(5 or 6)-Dimethyl-1,2-dithian-4-one (E or Z) | 4.59 | 3.98 | − | 8.6 | − | 12.26 | ms |
| 1266 | 3-Methyl-1,2,4-trithiane | 6.4 | − | 17.5 | − | 15.93 | − | MS |
| 1413 | 3,6-Dimethyl-1,2,4,5-tetrathiacyclohexane (E or Z) | − | − | 13.71 | 0.15 | 40.51 | − | MS |
| 1428 | 3,6-Dimethyl-1,2,4,5-tetrathiacyclohexane (E or Z) | − | − | 2.38 | − | 2.44 | − | MS |
| 1471 | 4,6-Dimethyl-1,2,3,5-tetrathiacyclohexane | 0.7 | − | 2.05 | − | 4.44 | − | MS |
| | *Disulfide* | | | | | | | |
| 1509 | 3-[(2-Methyl-3-furyl)dithio]-2-butanone | 36.62 | 3.19 | 8.57 | 11.92 | 4.85 | 3.67 | LRI |
| 1542 | bis(2-Methyl-3-furyl)disulfide | 3.94 | 7.87 | 1.23 | 10.39 | − | 1.4 | MS + LRI |
| 1578 | 3-[(2-Methyl-3-furyl)dithio]-2-pentanone | 0.82 | 2.25 | − | 3.97 | − | tr | LRI |
| 1592 | 2-[(2-Methyl-3-furyl)dithio]-3-pentanone | 2.76 | 1.46 | − | 2.32 | − | 1.8 | LRI |
| 1596 | 3-(2-Furfuryldithio)-2-butanone | 1.59 | 0.83 | 2.28 | 1.42 | 1.73 | tr | LRI |
| 1643 | 2-Methyl-3-(2-furfuryldithio)-furan | 0.48 | 0.85 | − | 1.49 | − | tr | LRI |
| 1668 | 3-(2-Furfuryldithio)-2-pentanone | − | 0.33 | 1.59 | 0.34 | 1.73 | 0.91 | LRI |
| 1702 | 2,3-Dihydro-5-methyl-4-[(2-methyl-3-furyl)dithio]furan | 3.2 | 6.84 | − | 8.82 | − | 1.14 | ms |
| 1745 | 2-Methyl-3-[(2-methyl-3-thienyl)dithio]furan | 3.93 | 2.73 | − | 4.17 | − | 0.37 | LRI |
| | *Thiazole* | | | | | | | |
| 933 | 4,5-Dimethylthiazole | 1.16 | − | 1.02 | − | − | − | ms + LRI |
| 997 | 2,4,5-Trimethylthiazole | + | + | 5.70 | + | 3.55 | + | ms + LRI |
| 1022 | 2-Acetylthiazole | 11.13 | 1.32 | 20.32 | 3.00 | 28.80 | 3.14 | MS + LRI |
| 1077 | 2,4-Dimethyl-5-ethylthiazole | 1.94 | − | 0.72 | 0.91 | 1.19 | 1.06 | ms + LRI |
| 1113 | 2-Acetyl-4-methylthiazole | − | − | 3.98 | − | 5.38 | − | ms |
| 1241 | 2,5-Diethyl-4-methylthiazole | − | − | 3.14 | − | 5.10 | 1.11 | ms |
| | *Pyrazine* | | | | | | | |
| 825 | Methylpyrazine | 5.75 | − | 21.95 | − | 37.97 | − | MS + LRI |
| 914 | 2,6-Dimethylpyrazine | − | − | 7.36 | − | 7.80 | − | MS + LRI |
| 918 | Ethylpyrazine | 2.87 | − | 3.24 | 0.91 | 6.51 | 1.51 | ms + LRI |
| 919 | 2,3-Dimethylpyrazine | 2.56 | − | 3.37 | 0.86 | 2.86 | 1.31 | MS + LRI |
| 999 | 2-Ethyl-6-methylpyrazine | 1.17 | − | 4.64 | − | 7.28 | − | ms + LRI |
| 1003 | Trimethylpyrazine | − | − | 1.89 | − | 1.67 | − | MS + LRI |
| | *Miscellaneous* | | | | | | | |
| 818 | 3-Mercapto-2-butanone | 6.88 | 4.65 | − | 4.13 | − | 0.66 | MS + LRI |
| 832 | Furfural | 3.19 | 14.29 | 1.74 | 7.17 | − | 3.80 | MS + LRI |
| 850 | 2-Furanmethanol | 2.29 | 0.57 | 1.14 | 1.35 | − | 1.60 | MS + LRI |
| 904 | 3-Mercapto-2-pentanone | 4.75 | 10.22 | − | 3.31 | − | 0.47 | MS + LRI |
| 909 | 2-Mercapto-3-pentanone | 6.32 | 3.02 | − | − | − | − | ms + LRI |
| 952 | 1-(2-Furyl)-2-propanone | 2.75 | 1.28 | − | 1.26 | − | 1.41 | ms |
| 986 | 1-Mercapto-3-pentanone | 2.6 | 1.92 | − | 1.50 | − | − | ms |
| 1009 | 1-(2-Furanyl)-1-propanone | 3.62 | − | 1.62 | − | − | 0.81 | ms |

Footnotes: ᵃ "Approximate quantities in headspace (ng mL⁻¹ of mixture) given as means of independent
experiments; tr, trace (<0.3 ng mL⁻¹ of mixture); −, below detection limit (~0.07 ng mL⁻¹ of mixture); +,
present in small amounts and quantification confounded by adjacent peak." ᵇ "Each model system consisted of
5 mmol ribose and 5 mmol cysteine in 100 mL 0.2 of mol L⁻¹ phosphate buffer." ᵈ MS + LRI = authentic
compound in the authors' lab; ms + LRI = library spectrum and LRI near literature; MS = literature spectrum;
ms = library spectrum only. Superscript literature numbers on the ID column condensed.

---

## 4. DERIVED [Z]

**Overlap with `xu2008_extraction.md` (FRI 41:730, 1:1 arm, initial pH 5.6):** MFT 16.62 / 34.42, FFT
22.09 / 25.13, 3-thiophenethiol 25.55 / 36.94, 2-methyl-3-thiophenethiol 25.37 / 14.22, 1325-dihydro-
thienothiophene 441.52 / 62.25, 2,3-dihydro-6-methylthieno[2,3-c]furan 56.08 / 13.13, methylpyrazine 5.75 —
**identical**; a few cells differ in the last digit (3-methyl-1,2-dithian-4-one 8.97 vs 8.96; 2-acetylthiazole
11.13 vs 11.12; furfural SC 14.29 vs 14.28; trithiolane SC "−" vs "tr"). **The acidic columns here are the
same experiment as xu2008's 1:1 columns, not a replicate.** New here: the pH 7.0 and 8.0 arms.

**pH ratios, conventional (final pH):** MFT 5.5 → 6.2: **6.3× fall**, → 6.8: > 237× (to < 0.07);
FFT 5.5 → 6.2: > 315×. 3-Mercapto-2-butanone/-pentanone also vanish at pH ≥ 6.2. Pyrazines rise
(methylpyrazine 5.75 → 21.95 → 37.97).

**Between-paper reproducibility of the nominal control pot** (50 + 50 mM, pyrophosphate pH 5.6, 140 °C,
1 h, same apparatus): this paper / xu2008 MFT 16.62, FFT 22.09 vs `xu2010` control MFT 7.2 ± 2.0, FFT
7.3 ± 1.9 ng/mL — **×2.3 (MFT) and ×3.0 (FFT) between two runs of the same protocol in the same lab.**

**Constant-pH comparison:** Meynier & Mottram 1995 (`meynier1995_extraction.md`, 140 °C, 1 h, pyrophosphate,
constant pH) MFT 15.2 → 1.0 µg/10 mg ribose from pH 4.5 to 6.0 (15×) and FFT 8.5 → 1.5 (5.7×); here MFT
falls 6.3× and FFT > 300× between final pH 5.5 and 6.2. Direction agrees; FFT magnitude much steeper here
(headspace SPME of a thiol near its detection limit).

---

## 5. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

- **pH 5.6 MFT/FFT row: already held** via `xu2008` (do not double-count).
- **New candidate HOLD-OUT shape rows:** MFT and FFT at **initial pH 7.0 and 8.0 (final 6.2, 6.8)**,
  140 °C / 1 h, 50 + 50 mM, 0.2 M pyrophosphate — a **three-point pH series with measured final pH** in the
  regime the core already models (140 °C). Levels are **headspace-equivalent ng/mL with RF = 1, not
  aqueous concentrations**, so only within-study ratios (MFT pH 5.5 : 6.2 = 6.3; FFT > 300) are usable,
  and the ×2–3 between-run spread (above) bounds what a "level" means. Not benchmark quality for absolute
  levels.
- **Disulfide sink:** the MFT mixed disulfides (with 3-mercapto-2-butanone: 36.62 ng/mL, > 2× free MFT at pH
  5.6) carry most of the MFT-derived signal, as in `mottram1995_extraction.md`; supports a mixed-disulfide
  sink beyond homodimers. Semi-quantitative.
- **SC-CO2 arm:** out of scope (pressure, CO2 acidification, CO2-stripping during depressurisation — the
  volatiles of the SC arm include three absorber traps the control does not have [!]).

## 6. FLAGS

1. **Same data as `xu2008`** for the pH 5.6 arms (§4).
2. **"Phosphate" buffer is pyrophosphate** by the Methods' own recipe.
3. **Year**: printed 2008 (vol. 88), online 2007; stem `xu2007`.
4. **Control arm gas/pressure not stated**; the SC arm sums reaction mixture + three absorbers, the control
   (presumably) only the mixture — the SC/control comparison is not like-for-like.
5. RF = 1 for all compounds; "not … absolute concentrations in the aqueous solutions" (authors).
6. MFT identified "ms + LRI" here but "MS + LRI" (authentic) in `xu2008`; FFT authentic in both.
7. Fig. 2 browning data not digitised.

**Relation to existing dossiers:** `xu2008_extraction.md` (ratio axis; shares the pH 5.6 data),
`xu2010_extraction.md` (pressure/atmosphere axis, same apparatus; ×2–3 different control), `meynier1995_extraction.md`
(constant-pH MFT/FFT at 140 °C), `mottram2002_extraction.md` (140 °C / 30 min pots), `mottram1995_extraction.md`
(mixed disulfides), `sun2022_extraction.md` (citing review).
