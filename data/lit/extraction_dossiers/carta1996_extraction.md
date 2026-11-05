# Carta & Tola 1996 — EXTRACTION
DOI: 10.1021/je9501853 (the PDF prints the manuscript code "JE9501853" on p. 417; Crossref 2026-11-05: "Solubilities of L-Cystine, L-Tyrosine, L-Leucine, and Glycine in Aqueous Solutions at Various pHs and NaCl Concentrations", Carta, Tola, J. Chem. Eng. Data 41(3), 414-417, 1996. **Matches the brief and the PDF.**)
### Equilibrium solubility of L-CYSTINE (with L-tyrosine, L-leucine, glycine) at **25 °C only**, pH 0.7-13 (HCl / NaOH), in water and in 1 and 3 M NaCl. **L-cystine: 0.69-0.71 mmol/L over initial pH 3.0-8.6 in water** (0.168 g/L), 0.83-0.87 mmol/L in 1 M NaCl, 1.09-1.10 mmol/L in 3 M NaCl; steep rise below pH 2 and above pH 10. Zwitterion solubility constant 6.92 × 10⁻⁴ M; salting-in constant +0.068 L/mol. **No temperature series, no phosphate.**

**Read 2026-11-05 from `data/articles/carta1996.pdf`** (4 pp., J. Chem. Eng. Data pp. 414-417; born-digital ACS PDF).
Tables 1-6 from `pdftotext -layout`; the text layer is clean (columns checked against the printed layout; no figure
values are used).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Dipartimento di Ingegneria Chimica e Materiali, Università di Cagliari (Carta; chemical
engineering, keratin-hydrolysate recovery). Not in `lab_ids.yml`; proposed id `cagliari_carta`. Independent of every
Maillard group in the corpus.

## 1. What this paper answers for the model

The question (case M1c sec. 4.3; the Hou 2026 dossier's flag): does cystine precipitate in the pots, and what
ceiling does solubility put on the cystine member P3 as it is measured?

| model need | what Carta 1996 gives | tag |
|---|---|---|
| L-cystine solubility, water, near-neutral | **0.69-0.71 mmol/L** (initial pH 3.04-8.61), 25 °C | [M] |
| salt effect | 1 M NaCl: 0.83-0.87 mmol/L (pH 4.2-9.3); 3 M NaCl: 1.09-1.10 mmol/L (pH 4.0-8.6); log(s/s0) = 0.068 C | [M] |
| pH effect | rises to 17.0 mM at pH 1.00 and 28 mM at pH 13 (water) | [M] |
| temperature effect | **none measured** (25 °C only) | — |
| phosphate / buffer effect | **none measured** | — |
| dissolution kinetics | cystine and tyrosine "very long dissolution times (about 24 h)"; 48 h equilibration | [Q] |

## 2. Conditions

| item | value |
|---|---|
| Solid | L-cystine, Aldrich Chemie 99 %, anhydrous (checked by mass after 48 h contact with water) |
| Solvent | bidistilled water + HCl (1 M stock) or NaOH (0.1 M stock) to a set pH; NaCl 0, 1, 3 mol/dm³ |
| Equilibration | excess solid (≥ 2x), magnetic stirring, **48 h**, **25.0 ± 0.1 °C**; approach from 20 °C and from 30 °C gave no significant difference |
| Sampling | ~10 cm³ by syringe through a thermostated filter |
| Analysis | UV at 215 nm (cystine), blank matched in acid / base; 1 : 2 dilution; ε 6.6778 dm³ g⁻¹ (r 0.998) |
| pH | Metrohm 691 (0.01 units), **of the initial solution, not of the saturated one** [Q]; near neutrality the two differ strongly (authors) |
| Replicates | three flasks × three samples; maximum deviation ± 1 × 10⁻⁵ mol/dm³ for L-cystine |

## 3. Tables re-typed (L-cystine columns only)

### Tables 1-3 (p. 415). Solubility s (mol/dm³), 25 °C [M]

| water: pH | s | 1 M NaCl: pH | s | 3 M NaCl: pH | s |
|---|---|---|---|---|---|
| 1.00 | 0.017 01 | 1.00 | 0.023 65 | 0.67 | 0.039 67 |
| 1.41 | 0.004 09 | 1.41 | 0.006 17 | 0.77 | 0.025 00 |
| 1.82 | 0.001 80 | 1.93 | 0.002 12 | 1.18 | 0.007 83 |
| 2.29 | 0.000 90 | 2.47 | 0.001 16 | 1.68 | 0.003 08 |
| 3.04 | 0.000 70 | 2.96 | 0.000 95 | 2.62 | 0.001 28 |
| 3.42 | 0.000 70 | 4.16 | 0.000 85 | 2.91 | 0.001 13 |
| 4.02 | 0.000 70 | 5.50 | 0.000 84 | 3.97 | 0.001 10 |
| 5.08 | 0.000 71 | 5.78 | 0.000 84 | 5.00 | 0.001 10 |
| 7.00 | 0.000 69 | 7.00 | 0.000 85 | 7.00 | 0.001 10 |
| 8.61 | 0.000 70 | 8.79 | 0.000 83 | 8.61 | 0.001 09 |
| 10.71 | 0.001 13 | 9.28 | 0.000 87 | 10.61 | 0.001 61 |
| 11.29 | 0.002 10 | 11.11 | 0.001 34 | 11.50 | 0.002 77 |
| 11.75 | 0.004 25 | 11.89 | 0.003 30 | 12.02 | 0.003 96 |
| 12.13 | 0.008 79 | 12.52 | 0.007 83 | 12.69 | 0.008 08 |
| 12.65 | 0.026 21 | 12.91 | 0.019 67 | 12.91 | 0.022 00 |
| 12.96 | 0.027 83 | 13.00 | 0.033 92 | 13.00 | 0.042 12 |
| 13.00 | 0.028 25 | | | | |

### Tables 4-6 (pp. 416-417) [M]

- K8 (zwitterion solubility constant, from the isoelectric data, γz = 1): **6.92 × 10⁻⁴ mol/dm³**.
- Dissociation constants used in the model (taken from Greenstein & Winitz 1961 [C]): K4 1.62 × 10⁻¹, K5 1.55 × 10⁻²,
  K6 5.75 × 10⁻⁹, K7 3.24 × 10⁻⁹ mol/dm³ (as printed; see Flag 3).
- Salting constant Kso (log(s/s0) = Kso·C_NaCl): **+0.068 dm³/mol** (salting in; tyrosine −0.044, leucine −0.112,
  glycine −0.002).
- Text: cystine solubility rises "by about 50 %" from 0 to 3 M NaCl.

## 4. Arithmetic (mine) [Z]

**1. Units.** 0.70 mmol/L × 240.3 g/mol = **0.168 g/L** at 25 °C in water; 0.85 mmol/L (1 M NaCl) = 0.204 g/L.
As cysteine equivalents: 1.40 mmol/L Cys bound as dissolved cystine at saturation (water, 25 °C).

**2. Salt check.** 3 M: log(1.10/0.70)/3 = 0.065; 1 M: log(0.845/0.70) = 0.082; the paper's regression 0.068. A 0.2 M
phosphate pot (ionic strength ~0.3-0.5 M at pH 5.5-6.5) would, if phosphate salted in like NaCl (untested), raise
the 25 °C solubility by ~5-8 % to ~0.74-0.76 mmol/L.

**3. The pots at 25 °C (cooled samples), against P3.**
- Zhang 2025b (Cys 40 mM): the case's P3 band holds ≤ ~1.4 mM Cys as cystine by 50 min = ≤ 0.7 mM cystine:
  **at, not above, the 25 °C saturation** (0.70-0.76 mM). Precipitation in the cooled sample is possible only at
  the band's ceiling, and would be slow (24-48 h dissolution kinetics; nucleation from supersaturation is not
  measured here).
- Precipitated cystine would remove cysteine equivalents only. It cannot produce Zhang's EQUAL Cys and Xyl deficit.
- Cao 2017 / Zhao 2017b (Cys 0.4 M): their cysteine-loss ceilings (≤ 2-4 × 10⁻⁴ min⁻¹) allow up to ~20 mM Cys lost
  over 5 h; if that went to cystine (≤ 10 mM) it would exceed 25 °C solubility ~14x. Their pools over-close, so
  nothing says it did.
- Hou 2026 storage arm (cystine 0.20-0.51 mM, read at room temperature): below saturation throughout. Hou's model pot
  (5 mM cystine charged): **7.1x** the 25 °C solubility (the Hou dossier said ~0.5 mM and "about tenfold").

## 5. What it can and cannot test in the model

- **Can:** set the room-temperature solubility ceiling (0.70 mM in water; ~0.75 mM at the pots' salt level if
  phosphate behaves like NaCl) that any measured dissolved cystine in a cooled sample cannot exceed at equilibrium;
  say that between pH 3 and 8.6 the ceiling is flat (zwitterion plateau).
- **Cannot:** say whether cystine precipitates AT process temperature (no T data; solubility of amino acids rises
  with T in general, not shown here), how fast it nucleates, or anything about phosphate.

## 6. Conflicts with existing records

| record | what it says | what Carta says | verdict |
|---|---|---|---|
| `hou2026_extraction.md` sec. on the cystine pot ("~0.11 g/L, ~0.5 mM, my value, not in the paper"; 5 mM "about tenfold" over) | 0.5 mM | **0.70 mM (0.168 g/L)** at 25 °C in water | **conflict in the number, not the conclusion**: Hou's cystine pot is 7x over saturation at room temperature, not 10x. The Hou dossier's value should be corrected by a later edit with this dossier as the source (not done here) |
| Dalton & Schmidt 1933 (cited in Fig. 1) | literature cystine solubility | agrees (Fig. 1, authors' statement) | [Q] |

## 7. Flags

1. **[!] pH is the INITIAL solution pH**, not the saturated solution's; near neutrality the authors say the two
   differ strongly. The flat 0.70 mM plateau is the isoelectric (zwitterion) solubility; the exact pH of each row
   between 3 and 8.6 is not known.
2. **[!] 25 °C only; no buffer salts.**
3. **[!] Table 5's constants as printed** (K4 1.62 × 10⁻¹, K5 1.55 × 10⁻², K6 5.75 × 10⁻⁹, K7 3.24 × 10⁻⁹) do not
   order as a tetraprotic acid's should for the labels given (K7, defined for A²⁺ → A⁺, is the smallest); they are
   literature values used in the authors' fit, not measured here. Not used.
4. The UV assay at 215 nm reads total dissolved cystine; it cannot see cysteine formed by any reduction.

## M1c role (frozen bd641ae8)

**Role fixed by the case (sec. 4.3):** band source for P3's ceiling, and for whether cystine precipitates on cooling
(an operator question).

**Band it can supply:**
- **Coordinate: P3 (Cys → ½ cystine) as seen by an assay of the cooled sample.** Dissolved cystine ≤ **0.70 mmol/L**
  (≤ 1.40 mmol/L Cys equivalents) at 25 °C in water, flat over pH 3-8.6; ~0.74-0.76 mmol/L at a 0.2 M phosphate pot's
  ionic strength if phosphate salts in like NaCl (Kso 0.068 L/mol; untested).
- **What it means for the case (arithmetic, no engine):** the case's P3 band (≤ ~1.4 mM Cys in 50 min at 120 °C,
  i.e. ≤ 0.7 mM cystine) sits exactly at the room-temperature saturation. So P3 as banded does not precipitate on
  cooling except at its ceiling, and if it did, it would remove cysteine only: precipitation cannot be Zhang's 1 : 1
  missing pool. It does not widen or narrow the P3 band (which is a rate ceiling from Hou 2026); it adds an operator
  ceiling on measured dissolved cystine.
- **Transfer:** the compound itself; 25 °C only, so it applies to the cooled sample, not the hot pot.
- **Not decisive.** It measures nothing the case would call decisive.
