# Cao et al. 2017 (FRI) — EXTRACTION (L-cysteine + D-xylose 0.4 mol/L each, ± glycine at 1-24× Cys, 0.2 M phosphate pH 5.5, 90 °C, 0.5-5 h; free Cys, free Xyl, TTCA, Cys-Amadori, Gly-Amadori, pH, A294/A420 all TIME-RESOLVED; volatiles at 3 h as relative peak areas)

### THE ONLY POT IN THE CORPUS THAT MEASURES TTCA *FORMATION* FROM FREE CYSTEINE + XYLOSE TOGETHER WITH THE FREE CYSTEINE AND FREE XYLOSE IT EQUILIBRATES WITH — six time points, four species, one pot. It also shows that at 90 °C and pH 4-5 TTCA's main fate is isomerisation to the cysteine Amadori compound, a species the engine does not carry.

**Source on disk:** `data/articles/cao2017.pdf` (38 pp., Elsevier **accepted manuscript**, not the typeset
version; FRIN 6740). Read from the `pdftotext -layout` text layer for prose and methods; **Table 1 was
re-typed from 200 dpi page renders** (the text layer interleaves the twelve columns and is unusable);
**Figs. 3, 4 and 5 were digitised from 220-400 dpi page renders** (scratchpad
`img/cao5a.png`, `cao5b.png`, `cao5c.png`, `cao4ab-35.png`, `cao4b-35.png`, `cao_fig-34.png`). Every
digitised number is tagged **[D]** with its read error. Supporting Information (Tables S1, S2) is
**not on disk**.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of glycine on reaction of cysteine-xylose: Insights on initial Maillard stage intermediates to develop meat flavor" |
| Authors | Changchun Cao, **Jianchun Xie** (corresponding, xjchun@th.btbu.edu.cn), Li Hou, Jian Zhao, Feng Chen, Qunfei Xiao, Mengyao Zhao, Mengdie Fan |
| Affiliations | Beijing Advanced Innovation Center for Food Nutrition and Human Health / Beijing Laboratory for Food Quality and Safety, Beijing Technology and Business University (BTBU); Feng Chen at Clemson University |
| Venue | Food Research International (2017) **99**, 444-453 (page range as cited by Zhao 2019; the accepted manuscript carries no pagination). Received 10 Feb 2017, revised 2 Jun 2017, accepted 2 Jun 2017 |
| DOI as printed | **10.1016/j.foodres.2017.06.012** ✔ matches the brief |
| Funding | NSFC 31371838, 31671895; National Key R&D Program 2017YFD0400100 |
| Same lab (Xie group, BTBU) | `hou2017_extraction.md` (fed Gly-Amadori / TTCA / Cys-Amadori, 140 °C, same buffer and vessel); Zhao 2019 (`zhao2019_extraction.md`, the ¹³C follow-up); Gong 2016 (in Chinese, not on disk) |
| Other-lab companions | Zhai 2019 / 2020 / 2021 / 2023 (Jiangnan, Ho), whose "ARP" is this paper's **Cys-Amadori** |

## 1. Why it matters

**(a) It is the formation half of the TTCA story, which the engine does not have.** `sulfur.py`
carries TTCA only as a **fed reactant** (`r_ttca_cys`: TTCA → Cys + PENT; `r_ttca_deg`: TTCA → DPO +
fragments). There is **no Cys + PENT → TTCA step**. This paper shows that in a free cysteine +
xylose pot the two condense to TTCA *fast and almost completely*: from 0.40 mol/L each, **TTCA
reaches 0.356 mol/L (89 % of the cysteine) by the first sample at 0.5 h at 90 °C [D]**, while free
cysteine falls to 0.076 mol/L and free xylose to 0.080 mol/L [D]. Zhai 2019 (`zhai2019_extraction.md`)
shows the same thing at pH 7.4 and adds that ~14.5 mmol/L TTCA is present **before heating**. Any
engine pot that charges free cysteine with a pentose therefore treats as free a cysteine pool that
is, in reality, mostly sequestered in the thiazolidine.

**(b) It gives the first apparent equilibrium constant for TTCA ⇌ Cys + Xyl in the corpus.**
`zhai2023jafc_extraction.md` appendix A.4(a) concluded that "the K for this one does not exist in the
corpus". From this paper's own digitised Figs. 4 and 5 (blank pot), the mass-action ratio
K_app = [TTCA] / ([Cys][Xyl]) is **58.5-64.6 L/mol at every one of the six sample times from 0.5 to
5 h** (arithmetic 1), while TTCA itself falls by a third and the pH drifts from 4.7 to 4.15. A
quotient that stays constant while all three concentrations move is the signature of a
**pre-equilibrium that is established within 30 min and then tracks the slow drain**. That is
exactly the missing reverse of `r_ttca_cys` that Zhai 2021 and Zhai 2023b inferred by sign only.

**(c) It shows TTCA's dominant fate at 90 °C / pH 4-5 is isomerisation, not ring-opening to free
Cys nor direct degradation to 1-DX.** In the blank pot TTCA falls 0.356 → 0.236 mol/L between 0.5
and 5 h while the Cys-Amadori rises 0.006 → 0.133 mol/L and free Cys stays at 0.06-0.08 mol/L [D]:
the cysteine-bearing pool (Cys + TTCA + Cys-Amadori) is conserved to within the analytical bias
(arithmetic 2). The authors' mechanism (after de Roos 1992, de Roos 2005) is TTCA ⇌ Cys-Amadori,
"acid-catalyzed by other amino acids", and the Cys-Amadori — not TTCA — is what degrades to the
meaty volatiles ("the rate-determining step to generate meat flavor should lie in the conversion
process from TTCA to Cys-Amadori but not the degradation of the latter"). **The engine has no
Cys-Amadori species**; its `ARP` is the *alanine* Amadori compound of Zhou 2023
(`species_sulfur.py`: N-(1-deoxy-D-xylulos-1-yl)-alanine), a different molecule.

**(d) Glycine accelerates the TTCA drain** (to ~0.03 mol/L by 2 h at 16-24× glycine [D]), which the
paper attributes to glycine catalysing TTCA → Cys-Amadori plus competition for xylose. Glycine
itself is not measured.

What this paper does NOT give: any fitted rate constant, any activation energy (one temperature
only), any volatile concentration (relative peak areas only), any H2S, any dicarbonyl, any
measurement at a second pH, or its SI.

### How it relates to the companion dossiers

| | **this paper** | `hou2017_extraction.md` (same lab) | `zhai2019_extraction.md` | `zhai2021_extraction.md` |
|---|---|---|---|---|
| start | **free Cys + Xyl** 0.4 M each (± Gly) | purified TTCA / Cys-Amadori / Gly-Amadori 30 mM | free Xyl 165 + Cys 82.6 mM | purified TTCA 10 mM (± Xyl) |
| medium | 0.2 M Na phosphate pH 5.5, 5 mL in 15 mL vial | **same buffer, same vessel** | unbuffered water, pH 7.4 (NaOH) | unbuffered, pH 5.5-8 |
| T, t | **90 °C**, 0.5-5 h | 140 °C, 5-30 min | 70-130 °C, 0-160 min | 100-140 °C, 20-120 min |
| TTCA measured | **yes, + Cys, Xyl, Cys-Amadori, pH** | no (UV, pH, volatiles) | yes, + Cys-Amadori | yes (+ Cys, Xyl, ARP) |
| Cys-Amadori / TTCA at end | **0.56** (5 h) | — | ≤ 0.06 (90 °C) | small |

Hou 2017 gives the **rank** Gly-Amadori > Cys-Amadori > TTCA in reactivity; this paper supplies the
**time-resolved concentrations** in the same buffer and vessel at a lower temperature.

## 2. Methods as they matter to a model

- **Charge.** "The amounts of cysteine and xylose used in the reaction models were both **2 mmol in
  5 mL** of phosphate buffer (**0.2 M, pH 5.5**)" → **0.40 mol/L each**, confirmed by the t = 0 points
  of Fig. 4 (both at 0.40 mol/L). Glycine "at **0, 1, 8, 12, 16, or 24 times** of the cysteine
  concentration (mol/L)" → 0 / 0.4 / 3.2 / 4.8 / 6.4 / **9.6 mol/L** (see Flag 1).
- **Vessel.** **15 mL pressure glass vial** (Beijing Synthware), Parallel Synthesis Poly-block 4
  (HEL), "sealed and heated at **90 °C** while stirring for **0.5, 1, 2, 3, 4, or 5 h**"; cooled "by a
  stream of tap water". **Headspace 10 mL of air** (5 mL liquid in 15 mL; atmosphere not stated, so
  air). Same vessel and fill as Hou 2017.
- **pH.** Set to 5.5 by the buffer; **measured after heating** with a PHSJ-5 meter (Fig. 4c) — one of
  few pots in the corpus with a measured pH trajectory.
- **Cys and Xyl.** HPLC-ELSD (Agilent 1100 + SEDEX 75), XBridge Amide 4.6 × 150 mm 3.5 µm, isocratic
  ACN : 10 mM ammonium formate pH 6 = 25 : 75, 1 mL/min, 25 °C, 5 µL; external standards.
  Calibrations: **Cys y = 38046 x − 311.57, R² 0.9995 (0.01-0.1 mol/L); Xyl y = 61719 x − 587.1,
  R² 0.9995 (0.05-0.175 mol/L)**; x in mol/L, y peak area. Note the Cys range tops out at 0.1 mol/L
  and the Xyl range starts at 0.05 mol/L (Flag 3).
- **TTCA, Cys-Amadori, Gly-Amadori.** LC-MS (LCQ-DECA XP MAX ion trap, ESI+), same column,
  gradient 95→75 % ACN 0-5 min, 75→65 % 5-25 min, 0.5 mL/min; **SIM m/z 254 [M+H]⁺ + 236
  [M+H−H₂O]⁺ for TTCA and Cys-Amadori (same masses, separated by retention only), m/z 208 for
  Gly-Amadori**. External standards of the synthesised compounds, 0.05-0.5 mol/L: **TTCA y = 72.398 x +
  0.613, R² 0.9902; Cys-Amadori y = 57.752 x − 0.132, R² 0.9985; Gly-Amadori y = 44.513 x + 0.293,
  R² 0.9928**.
- **Standards.** Gly-Amadori by Hao 2007; Cys-Amadori and TTCA by de Roos 2005 (TTCA isolated from a
  **50 °C / 1 h** reaction); Bio-Rad AG 50W-X4 (H⁺) cation exchange, eluted with 0.3 M ammonia;
  purities by HPLC-ELSD peak area **> 99 % (Gly-Amadori), > 98 % (Cys-Amadori), > 99 % (TTCA)**;
  identity by MS and ¹³C NMR (D₂O, 300 MHz), spectra not printed.
- **UV.** A294 and A420 (Shimadzu UV 2700), diluted and **back-converted to the undiluted solution**.
- **Volatiles (3 h samples only).** 1 mL reaction mixture + 4 mL phosphate buffer + 1 µL
  1,2-dichlorobenzene (0.1 mg/mL in MeOH) in a 15 mL vial; 50/35 µm CAR/PDMS/DVB fibre (as printed), 10 min
  equilibration + 30 min exposure at **30 °C**; Agilent 7890A/5975C, HP-5MS and DB-Wax; quantities
  "approximate", **relative peak areas** to the internal standard (Farmer & Mottram 1990).
- **Statistics.** "All the results were averages of three replicates. The relative standard
  deviations (RSD) were all less than 10%." Duncan's test, p < 0.05.

## 3. Tables re-typed

### Table 1 (verbatim). Volatile flavor compounds in 3 h reaction mixtures of cysteine-xylose with glycine from zero (blank) to 24× the cysteine concentration (mol/L), SPME/GC-MS on HP-5 and DB-Wax

Amounts are **relative peak areas to 1,2-dichlorobenzene (dimensionless), mean ± SD, n = 3**;
superscript letters = Duncan groups within a row for one column; "-" not detected. Column headers are
the glycine : cysteine multiple.

| Compound | RI HP-5 | RI DB-Wax | HP-5 Blank | HP-5 1 | HP-5 8 | HP-5 12 | HP-5 16 | HP-5 24 | Wax Blank | Wax 1 | Wax 8 | Wax 12 | Wax 16 | Wax 24 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2-methylfuran | 621 | 888 | 0.08±0.00ᶜ | 0.18±0.02ᵇᶜ | 0.58±0.07ᵇ | 1.57±0.11ᵃ | 1.47±0.05ᵃ | 1.23±0.07ᵃ | - | 0.11±0.01ᶜᵈ | 0.65±0.11ᶜᵈ | 2.98±0.07ᵃᵇ | 5.22±0.22ᵃ | 1.31±0.06ᵇᶜᵈ |
| pyrazine | 730 | 1205 | - | - | - | - | - | - | 0.01±0.00ᶜ | 0.02±0.00ᵇᶜ | 0.10±0.00ᵃ | 0.05±0.01ᵃᵇᶜ | - | - |
| thiazole | 736 | 1244 | - | - | - | - | - | - | - | 0.02±0.00ᵇ | 0.04±0.00ᵃ | 0.02±0.00ᵇ | 0.01±0.00ᶜ | 0.02±0.00ᵇ |
| 2-methylthiophene | 774 | 1088 | - | - | 0.17±0.01ᵇ | 0.32±0.00ᵃᵇ | 0.26±0.03ᵃᵇ | 0.37±0.00ᵃ | - | - | 0.47±0.06ᵇᶜ | 2.14±0.08ᵃ | 3.16±0.11ᵃ | 0.18±0.01ᶜ |
| dihydro-2-methyl-3(2H)-furanone | 806 | 1255 | - | - | - | - | - | - | - | - | 0.03±0.00ᶜ | 0.12±0.00ᵇ | 1.41±0.20ᵃ | 0.13±0.01ᵇ |
| 3-mercapto-2-butanone | 816 | 1291 | - | - | 0.03±0.00ᵃ | 0.06±0.00ᵃ | 0.07±0.00ᵃ | 0.04±0.00ᵃ | - | - | 0.06±0.00ᵇ | 0.18±0.01ᵃ | 0.22±0.00ᵃ | 0.19±0.06ᵃ |
| furfural | 833 | 1456 | - | - | - | - | - | - | 0.06±0.00ᵈ | 0.07±0.00ᶜᵈ | 0.13±0.00ᵇᶜᵈ | 0.28±0.01ᵃ | 0.26±0.01ᵃᵇ | - |
| 2-methylthiazole | 852 | 1234 | - | - | - | - | - | - | - | - | 0.02±0.00 | - | - | - |
| **2-methyl-3-furanthiol** | 866 | 1304 | 0.38±0.04ᵉ | 3.71±0.43ᵉ | 51.82±2.61ᵇ | 58.14±2.17ᵃ | 35.07±3.98ᶜ | 22.03±1.00ᵈ | 0.02±0.00ᵉ | 0.39±0.02ᵉ | 21.90±2.58ᵈ | 47.58±1.13ᵃ | 28.74±2.33ᵇᶜᵈ | 24.17±0.30ᶜᵈ |
| 3-mercapto-2-pentanone | 901 | 1347 | 0.01±0.00ᶜ | 0.09±0.01ᶜ | 1.43±0.12ᵇ | 3.47±0.37ᵃ | 1.86±0.71ᵇ | 1.25±0.05ᵇ | - | 0.06±0.00ᶜ | 5.51±0.12ᵃ | 9.64±0.01ᵃ | 8.09±0.52ᵃ | 7.72±0.68ᵃ |
| **2-furfurylthiol** | 910 | 1428 | 0.02±0.00ᶜ | 0.04±0.00ᶜ | 0.18±0.02ᵇ | 0.66±0.00ᵃ | 0.26±0.03ᵇ | 0.24±0.00ᵇ | 0.03±0.00ᵉ | 0.03±0.00ᵉ | 0.55±0.02ᵇᶜ | 0.80±0.02ᵃ | 0.25±0.00ᵈᶜ | 0.36±0.03ᶜ |
| 2-methyl-3-(methylthio)furan | 947 | - | - | 0.02±0.00ᵃ | - | 0.01±0.00ᵃ | - | 0.02±0.00ᵃ | - | - | - | - | - | - |
| 2-thiophenethiol | 965 | 1530 | - | 0.01±0.00ᶜ | 0.13±0.01ᵃ | 0.12±0.00ᵃ | 0.04±0.00ᵇ | 0.04±0.00ᵇ | - | - | - | - | - | - |
| dihydro-2-methyl-3(2H)-thiophenone | 983 | 1512 | - | 0.01±0.00ᵇ | 0.02±0.00ᵇ | 0.04±0.00ᵇ | 0.04±0.00ᵇ | 0.09±0.00ᵃ | - | - | 0.11±0.00ᵈ | 0.19±0.00ᵇᶜ | 0.24±0.00ᵃ | 0.18±0.05ᵇᶜᵈ |
| 2-acetylthiazole | 1015 | 1631 | - | - | - | - | - | - | - | - | 0.05±0.00ᵃ | 0.12±0.00ᵃ | - | - |
| 2-[(methylthio)methyl]-furan | 1017 | 1542 | - | - | 0.04±0.00ᵇ | 0.07±0.01ᵃ | 0.09±0.01ᵃ | 0.02±0.00ᵇ | - | - | 0.06±0.00ᵈ | 0.17±0.00ᵇᶜ | 0.24±0.00ᵃ | 0.15±0.01ᵇᶜ |
| 1-(3-thienyl)-ethanone | 1044 | 1759 | - | - | - | - | - | - | - | - | 0.04±0.01ᵃ | - | - | 0.02±0.00ᵃ |
| 3-(methylthio)thiophene | 1049 | 1635 | - | - | - | 0.05±0.00 | - | - | - | - | - | - | - | - |
| 2-formyl-2,3-dihydro-thiophene | 1078 | 1336 | - | - | - | - | - | 0.02±0.00 | - | 0.10±0.00ᶜ | 2.69±0.10ᵇ | 4.62±0.25ᵃ | 2.86±0.12ᵇ | 2.35±0.14ᵇ |
| 3-methyl-2-thiophenecarboxaldehyde | 1116 | 1795 | - | - | 0.05±0.00ᵃ | 0.01±0.00ᵇ | - | - | - | - | 0.04±0.00ᵇᶜ | 0.09±0.00ᵃ | 0.06±0.00ᵇᶜ | 0.02±0.00ᶜ |
| 3,3'-dithiobis-(2-methyl)-furan | 1521 | 2156 | 0.07±0.01ᶜ | 0.32±0.04ᶜ | 6.46±0.34ᵃ | 7.64±0.78ᵃ | 4.15±0.28ᵇ | 2.81±0.14ᵇ | 0.03±0.00ᵉᶠ | 0.16±0.02ᵉᶠ | 5.98±0.58ᵇᶜᵈ | 13.71±0.35ᵃ | 9.27±1.45ᵇᶜ | 4.95±0.68ᵈ |
| 2-[(2-methyl-3-furyl)dithio]-3-pentanone | 1571 | 2426 | - | - | - | - | - | 0.03±0.00 | - | - | 0.02±0.00 | - | - | - |
| **Total** | | | 0.56 | 4.38 | 60.91 | 72.16 | 43.31 | 28.19 | 0.15 | 0.96 | 38.45 | 82.69 | 60.03 | 41.75 |

Footnotes (verbatim in substance): ¹ identified by NIST 11, RI against published values and
injection of authentic chemicals, **except** 2-methylfuran, 2-formyl-2,3-dihydro-thiophene,
2-[(methylthio)methyl]-furan, 2-methyl-3-(methylthio)furan and 2-[(2-methyl-3-furyl)dithio]-3-pentanone
(tentative); ² linear RIs vs C6-C26 n-alkanes; ³ means ± SD (n = 3), relative peak areas to
1,2-dichlorobenzene; different letters in a row on one column = p < 0.05; "-" not detected.

The 2-methyl-3-furanthiol row is the column-by-column reading of the 200 dpi render; the
superscript on its DB-Wax 24× cell reads "cd".

### Digitised Fig. 5 — TTCA, Cys-Amadori, Gly-Amadori against time, 90 °C [D]

Axis 0-0.40 mol/L; read error **± 0.005 mol/L** (≈ ± 8 px at 400 dpi); error bars visible only on
some points and not re-read. t = 0 is 0 for all species in every arm.

**Fig. 5a TTCA (mol/L)**

| arm (Gly : Cys) | 0.5 h | 1 h | 2 h | 3 h | 4 h | 5 h |
|---|---|---|---|---|---|---|
| blank | 0.356 | 0.338 | 0.290 | 0.268 | 0.254 | 0.236 |
| ×1 | 0.315 | 0.304 | 0.283 | 0.261 | 0.244 | 0.227 |
| ×8 | 0.276 | 0.238 | 0.164 | 0.091 | 0.055 | 0.044 |
| ×12 | 0.234 | 0.195 | 0.073 | 0.051 | 0.038 | 0.033 |
| ×16 | 0.193 | 0.126 | 0.047 | 0.034 | 0.033 | 0.023 |
| ×24 | 0.158 | 0.081 | 0.029 | 0.020 | 0.016 | 0.013 |

The blank 0.5 h point carries the largest error bar in the figure (≈ 0.32-0.39 mol/L).

**Fig. 5b Cys-Amadori (mol/L)**

| arm | 0.5 h | 1 h | 2 h | 3 h | 4 h | 5 h |
|---|---|---|---|---|---|---|
| blank | 0.006 | 0.063 | 0.090 | 0.106 | 0.109 | 0.133 |
| ×1 | 0.005 | 0.060 | 0.085 | 0.097 | 0.103 | 0.099 |
| ×8 | ~0.005 | 0.052 | 0.080 | 0.083 | 0.079 | 0.069 |
| ×12 | ~0.005 | 0.049 | 0.072 | 0.063 | 0.054 | 0.046 |
| ×16 | ~0.005 | 0.044 | 0.049 | 0.053 | 0.038 | 0.036 |
| ×24 | ~0.005 | 0.038 | 0.044 | 0.035 | 0.026 | 0.023 |

**Fig. 5c Gly-Amadori (mol/L)** (no blank series)

| arm | 0.5 h | 1 h | 2 h | 3 h | 4 h | 5 h |
|---|---|---|---|---|---|---|
| ×1 | 0.003 | 0.010 | 0.010 | 0.010 | 0.007 | 0.006 |
| ×8 | 0.012 | 0.058 | 0.067 | 0.065 | 0.043 | 0.028 |
| ×12 | 0.024 | 0.083 | 0.083 | 0.057 | 0.039 | 0.031 |
| ×16 | 0.031 | 0.153 | 0.049 | 0.016 | 0.008 | 0.004 |
| ×24 | 0.038 | 0.176 | 0.073 | 0.043 | 0.026 | 0.016 |

**Read ambiguity (Fig. 5c):** at 1 h there are **six** markers for five series — a second triangle
sits at ≈ 0.165 mol/L, overlapping the ×24 square. The ×8 assignment above follows the connected
line (0.012 → 0.058 → 0.067 …); the high triangle is unassigned. The text's "maximal value of
Gly-Amadori arrived at around 1 h" holds for ×16 and ×24, not for ×8 or ×12.

### Digitised Fig. 4 — free cysteine, free xylose, pH, 90 °C [D]

Read error **± 0.004 mol/L** (Cys, Xyl), **± 0.03 pH units**. Both sugars and Cys are 0.40 mol/L at
t = 0 in all arms; pH at t = 0 is 5.42-5.49.

| species | arm | 0.5 h | 1 h | 2 h | 3 h | 4 h | 5 h |
|---|---|---|---|---|---|---|---|
| **Cys** | blank | 0.076 | 0.068 | 0.066 | 0.064 | 0.062 | 0.056 |
| Cys | ×1 | 0.068 | 0.062 | 0.059 | 0.056 | 0.050 | 0.043 |
| Cys | ×8-×16 (bundle, not separable) | ~0.050 | ~0.045 | ~0.037 | ~0.030 | ~0.026 | ~0.024 |
| Cys | ×24 | 0.046 | 0.032 | 0.030 | 0.026 | 0.023 | 0.019 |
| **Xyl** | blank | 0.080 | 0.077 | 0.073 | 0.071 | 0.070 | 0.069 |
| Xyl | ×1 | 0.069 | 0.064 | 0.062 | 0.060 | 0.056 | 0.054 |
| Xyl | ×24 | 0.022 | 0.014 | ~0.010 | ~0.010 | ~0.010 | ~0.010 |
| **pH** | blank (5.49 at 0 h) | 4.67 | 4.56 | 4.49 | 4.33 | 4.26 | 4.15 |
| pH | ×24 (5.42 at 0 h) | 4.43 | 4.24 | 4.11 | 4.07 | 4.05 | 4.03 |

### Digitised Fig. 3 — A294 and A420 (undiluted basis) [D], coarse (± 10 / ± 5 units)

A420 at 5 h: **blank ≈ 78**, ×1 ≈ 125, ×8 ≈ 333, ×12 ≈ 370, ×16 ≈ 403, ×24 ≈ 450. A294 at 5 h: blank
≈ 100, ×1 ≈ 270, ×8 ≈ 1560, ×12 ≈ 1770, ×16 ≈ 1950, ×24 ≈ 2010. The blank A420 is ≈ 0 to 3 h and
rises only after (≈ 45 at 4 h). See Flag 4 on the "ten times" claim.

### Arithmetic on the digitised values (all mine)

**1. Apparent equilibrium quotient of TTCA ⇌ Cys + Xyl, blank pot.**
K_app = [TTCA]/([Cys][Xyl]) in L/mol: **58.6 / 64.6 / 60.2 / 59.0 / 58.5 / 61.1** at 0.5 / 1 / 2 / 3 / 4 /
5 h. Mean **60 L/mol, spread ± 4 %**; read-error propagation alone gives ± 7 %, and the pool
over-closure (arithmetic 2) implies a systematic bias of up to ~15 %, so carry **K_app ≈ 60 L/mol
(plausible range 45-75)** at **90 °C, pH 4.2-4.7, ionic strength ~0.4 M phosphate, 0.4 M reactants**.
The ×1-glycine arm gives 67 → 98 L/mol, drifting upward with time — consistent with glycine
draining free Cys and Xyl through other channels, so **use the blank arm only**.
The stoichiometric K_app lumps every protonation state of Cys and TTCA; it is pH-dependent and
**cannot be carried to pH 7 without a speciation model** (TTCA ring N pKa 6.24 per the engine;
cysteine amine pKa ~8.3 at 25 °C). Evidence class **derived_from_digitised**.

**2. Cysteine-pool closure, blank.** Cys + TTCA + Cys-Amadori = 0.438 / 0.469 / 0.446 / 0.438 /
0.425 / 0.425 mol/L against 0.400 charged. The pool **over-closes by 6-17 %**, i.e. the three
calibrations are biased high together (TTCA R² is the worst, 0.9902), and then declines by
**~0.044 mol/L (≈ 10 %) between 1 and 5 h**. Read as an upper bound on *all* cysteine-destroying
channels (H2S release, cystine, Strecker, degradation of TTCA and Cys-Amadori to dicarbonyls) at
90 °C, pH 4.2-4.6: pooled first-order loss **≤ ~4 × 10⁻⁴ min⁻¹** (0.1/240 min). For scale, Kang
2026's free-cysteine consumption at 100 °C, pH 7, 10 mM is ~1.5 × 10⁻³ min⁻¹ (16.2 % in 120 min).

**3. TTCA → Cys-Amadori, blank.** Net Cys-Amadori gain 0.070 mol/L over 1-5 h against mean TTCA
≈ 0.28 mol/L → net first-order **~1.0 × 10⁻³ min⁻¹** (forward minus reverse minus Cys-Amadori loss,
so a **floor** on the forward isomerisation). Over 0.5-1 h the gain is 0.057 mol/L in 30 min →
~5.5 × 10⁻³ min⁻¹; the early burst may be Cys-Amadori forming directly from the open Schiff base
rather than through TTCA. Cys-Amadori / TTCA = **0.56 at 5 h and still rising** (the authors
report, data not shown, that it peaks at **7 h** in the blank and **5.2 h** in the ×1 arm), so the
TTCA ⇌ Cys-Amadori equilibrium, if there is one, lies at **≥ 0.56** at 90 °C and pH ~4.2.

**4. A floor on the condensation rate.** If TTCA is ≥ 90 % equilibrated by 30 min, the relaxation
rate k_f(2[Cys]_eq) + k_r ≥ ln 10/30 min = 0.077 min⁻¹, with [Cys]_eq ≈ 0.078 mol/L and K = 60 L/mol:
**k_f ≥ ~0.44 L mol⁻¹ min⁻¹ and k_r = k_f/K ≥ ~7 × 10⁻³ min⁻¹** at 90 °C, pH ~4.7. The ring-opening
is therefore *fast* compared with every net TTCA loss in the corpus (Zhai 2021's first-order
equivalents 0.0027-0.0081 min⁻¹ at 100-140 °C, pH 7) — the net loss is slow only because the
released Cys and Xyl recombine. **This is the quantitative form of the B16 diagnosis** recorded in
`zhai2023jafc_extraction.md` (appendix A): a one-way `r_ttca_cys` fitted where the reverse is active cannot be right.
Class **derived_assumption** (assumes the 0.5 h point is equilibrated; the 0 → 0.5 h interval is
unsampled).

**5. Glycine solubility.** 24 × 0.40 = **9.6 mol/L glycine = 720 g/L**; 8× is already 3.2 mol/L
(240 g/L). Glycine's solubility in water is ~250 g/L at 25 °C; the high-glycine arms cannot have
been single-phase solutions at charging (Flag 1).

## 4. Kinetic numbers the repository can use

Shared conditions: Cys 0.40 mol/L + Xyl 0.40 mol/L in 0.2 M sodium phosphate, initial pH 5.5
(measured 5.42-5.49), 5 mL in a sealed 15 mL glass vial under air, stirred, 90 °C, tap-water quench,
n = 3, RSD < 10 %.

| step / observable | quantity | value | unit | conditions | order as fitted | source | evidence class |
|---|---|---|---|---|---|---|---|
| TTCA formation from free Cys + Xyl | TTCA vs t, 6 arms | Fig. 5a table | mol/L | 90 °C, 0.5-5 h | **nothing fitted in the paper** | Fig. 5a [D] | **digitised time series** (± 0.005) |
| free Cys | Cys vs t | Fig. 4a table | mol/L | same | — | Fig. 4a [D] | digitised time series (± 0.004) |
| free Xyl | Xyl vs t | Fig. 4b table | mol/L | same | — | Fig. 4b [D] | digitised time series (± 0.004) |
| Cys-Amadori (xylose-cysteine Amadori) | vs t | Fig. 5b table | mol/L | same | — | Fig. 5b [D] | digitised time series — **no engine species** |
| Gly-Amadori | vs t | Fig. 5c table | mol/L | glycine arms | — | Fig. 5c [D] | digitised; one marker ambiguous |
| pH trajectory | pH vs t | 5.49 → 4.15 (blank) | — | same | — | Fig. 4c [D] | digitised (± 0.03) |
| **TTCA ⇌ Cys + Xyl** | K_app = [TTCA]/([Cys][Xyl]) | **≈ 60 (45-75)** | L/mol | 90 °C, pH 4.2-4.7, blank | equilibrium | arithmetic 1 | **derived_from_digitised** |
| TTCA ring-opening (reverse of condensation) | k_r floor | **≥ ~7 × 10⁻³** | min⁻¹ | 90 °C, pH ~4.7 | 1st order | arithmetic 4 | derived_assumption |
| TTCA → Cys-Amadori | net first-order floor | **~1.0 × 10⁻³** | min⁻¹ | 90 °C, pH 4.2-4.6, blank 1-5 h | 1st order | arithmetic 3 | derived_from_digitised |
| Cys-Amadori / TTCA | ratio at 5 h | **0.56**, rising | — | blank | — | Fig. 5 [D] | within_study_ratio |
| all cysteine-destroying channels pooled | first-order ceiling | **≤ ~4 × 10⁻⁴** | min⁻¹ | 90 °C, pH 4.2-4.6, blank | — | arithmetic 2 | derived_from_digitised (bound) |
| browning | A420 at 5 h, ×24 / blank | ~5.8 (text says ~10) | — | 90 °C | — | Fig. 3b [D] | within_study_ratio |
| volatiles at 3 h | relative peak areas | Table 1 | dimensionless | 90 °C, 3 h | — | Table 1 | **level_only** (relative areas, no calibration) |
| MFT optimum vs glycine | bell-shaped, max at 12× | Table 1 | — | 3 h | — | Table 1 | within_study_ratio |

## 5. What this can and cannot anchor in the model

**Can (with structural changes the engine does not yet have):**
1. **A TTCA formation step and its equilibrium.** The blank arm is a clean four-species,
   six-time-point record of Cys + Xyl ⇌ TTCA at one temperature, with measured pH. It is the
   natural **fit row for k_f/k_r of a reversible condensation** and fixes K_app ≈ 60 L/mol at
   90 °C / pH ~4.5. With Zhai 2019 (pH 7.4, 70-130 °C) it would bound the temperature dependence.
2. **The structural argument against a one-way `r_ttca_cys`.** Arithmetic 4 puts the ring-opening at
   ≥ 7 × 10⁻³ min⁻¹ at 90 °C while the net TTCA drain is ~10× slower. A one-way step fitted to a
   free-cysteine reading in a pot where recombination is active will overestimate net TTCA loss —
   which is B16's measured residual (~10× too fast on `k_ttca_cys`).
3. **A missing isomer.** The blank arm says ≥ 90 % of the TTCA lost over 0.5-5 h becomes
   Cys-Amadori at 90 °C / pH 4-5. A `CysARP` species with TTCA ⇌ CysARP and CysARP → DPO/TDP would
   be needed to score anything from this pot beyond TTCA itself.
4. **A cysteine-sink ceiling** (arithmetic 2) usable as a hold-out on the lumped `r_cys_thermal`
   and friends at 90 °C, pH 4-5.

**Cannot:**
- No volatile row: Table 1 is relative peak areas at one time; no calibration.
- No temperature dependence: one temperature.
- No pH-7 transfer of K without a speciation model.
- The glycine arms are not clean solutions (Flag 1) and glycine is not measured.
- Nothing on H2S, dicarbonyls, MFT/FFT kinetics or any thiol-consumption channel.

## 6. Flags

1. **The high-glycine arms are physically implausible as written.** 8×-24× glycine at 0.40 mol/L
   cysteine is 3.2-9.6 mol/L (240-720 g/L). These cannot be homogeneous solutions at charging; they
   are slurries that dissolve (partly) on heating, with water activity well below the blank. Treat
   every glycine arm as qualitative.
2. **TTCA and Cys-Amadori share m/z 254/236** and are separated only by retention in SIM LC-MS. A
   co-elution or response crossover would shift mass between them. The pool over-closes by
   6-17 % (arithmetic 2); the bias is systematic, not random.
3. **Calibration ranges vs. measured values.** Cys calibrated 0.01-0.1 mol/L, but the t = 0 value
   (0.40) is 4× above range (it is presumably the charge, not a measurement). Xyl calibrated
   0.05-0.175 mol/L, so the glycine-arm xylose values (0.010-0.022) are **below** the lowest
   standard. TTCA/Cys-Amadori/Gly-Amadori calibrated 0.05-0.5 mol/L; every Cys-Amadori and
   Gly-Amadori point below 0.05 mol/L, and every high-glycine TTCA point after 2 h, is an
   extrapolation.
4. **"About ten times" is not what Fig. 3b shows.** The text says A420 at 5 h for 24× glycine "was
   about ten times of that of the cysteine-xylose alone"; the figure reads ≈ 450 vs ≈ 78, i.e.
   **≈ 5.8×** [D].
5. **Fig. 5c has one marker too many at 1 h** (six markers, five series); see the read-ambiguity
   note.
6. **Accepted manuscript, not the version of record.** Pagination 99:444-453 is from the citing
   paper (Zhao 2019 ref list). Numbers could differ from the typeset article.
7. **"Cys-Amadori" ≠ the engine's `ARP`.** This paper's Cys-Amadori is
   N-(1-deoxy-D-xylulos-1-yl)-L-cysteine; the engine's `ARP` is the alanine Amadori compound. Do not
   map one onto the other.
8. **pH is buffered but not held.** 0.2 M phosphate at pH 5.5 against 0.4 M cysteine is a weak
   buffer (phosphate pKa₂ ≈ 6.8 sits outside the working range) and the pH drifts by 1.3 units in
   5 h. Any K or rate here is for a **drifting** pH.
9. **High concentration.** 0.4 mol/L reactants (40× the TTCA pots of Kang/Zhai). Activity
   corrections on K_app are not negligible; transfer to 10 mmol/L pots is order-of-magnitude.
10. **Registry gaps** (`data/keys/compounds.yml`): TTCA, Cys-Amadori, Gly-Amadori, cysteine, xylose,
    glycine are absent; `2_methyl_3_furanthiol`, `2_furfurylthiol`, `2_methylthiophene`, `furfural`
    are present. Paper not in `data/keys/papers.yml`.
11. **What to request from the authors:** the numeric data behind Figs. 3-5 (with SDs); the 7 h
    extension of the blank and ×1 arms ("data not shown"); glycine concentrations; Tables S1/S2
    (volatiles at the other time points).
