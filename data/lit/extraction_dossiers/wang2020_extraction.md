# Wang T., Zhen, Tan, Xie, Cheng & Zhao 2020 — EXTRACTION
### *"Characterization of initial reaction intermediates in heated model systems of glucose, glutathione, and aliphatic aldehydes"*
### Glutathione + glucose 60 + 60 mmol/L ± hexanal, (E)-2-heptenal or (E,E)-2,4-decadienal 140 mmol/L in 0.2 M phosphate pH 6.5, sealed tube, 140 °C, 0.16–4 h: time courses (Fig. 1h–j) of the GSH-Amadori compound, three Cys-Gly–sugar thiazolidines and ten aldehyde–GSH / aldehyde–Cys-Gly adducts, by HPLC-MS against GSH as the calibrant. Structures by MS/MS and NMR. Volatiles from re-heated purified adducts only as relative peak area (Table S1, not on disk). No free GSH time course, no H₂S, no MFT/FFT numbers. Read 2026-10-27.

**Source on disk:** `data/articles/wang2020.pdf` (publisher's PDF, 10 pp.). Text from `pdftotext -layout`;
Fig. 1(h–j) (p. 4) read from a 300 dpi render, panels cropped and enlarged. Supplementary material
(Fig. S1 scheme, Figs. S2–S6 NMR, Table S1 volatiles) **not on disk**.

**Provenance:** **[M]** printed number · **[Q]** qualitative statement · **[F]** read off a figure by me
(read error stated) · **[Z]** derived by me, never printed · **[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | Tianze Wang, Dawei Zhen, Jia Tan, Jianchun Xie\*, Jie Cheng, Jian Zhao |
| Title | "Characterization of initial reaction intermediates in heated model systems of glucose, glutathione, and aliphatic aldehydes" |
| Venue | *Food Chemistry* **305** (2020) **125482** |
| DOI | 10.1016/j.foodchem.2019.125482 (printed) |
| Affiliation | Beijing Technology and Business University; Beijing Lanjingzhongyu Scientific Development Co.; CAAS |
| Dates | received 16 Mar 2019; revised 4 Sep 2019; accepted 4 Sep 2019; online 5 Sep 2019 |

Follow-up to `zhao2019b_extraction.md` (same group, "reaction conditions were used according to our
former investigation" — i.e. the 140 °C GSH + glucose pot). Cited by `zhou2023c_extraction.md` (ref 20)
for "GSH degraded into Cys first, enhancing H₂S release". Related: `du2025_extraction.md` (GSH-Amadori
time course at 145 °C), `sun2019b_extraction.md` (GSH-xylose Amadori).

---

## 1. THE ONE-PARAGRAPH ANSWER

At 140 °C, pH 6.5 (buffer label), 60 mM GSH + 60 mM glucose ± 140 mM aldehyde, **every glutathione-derived
intermediate is transient on a 10–60 min scale**. The **GSH-Amadori compound (G4) peaks at the first
sample (~0.16 h) at ~0.175 mmol/L — 0.3 % of the GSH charge — and is essentially gone by 0.5 h**.
Glucose-side products that persist longer are **thiazolidines of Cys-Gly** (GSH minus its γ-glutamyl)
with intact glucose (G3, peak ~0.26 mM at 0.25 h) or C3/C4 sugar fragments (G1, G2, peak 0.14–0.28 mM at
0.5 h). Cys-Gly itself is a major peak already at 0.16 h in the GSH + glucose control chromatogram (Fig.
1a), i.e. **γ-glutamyl loss is fast**. With α,β-unsaturated aldehydes GSH is captured far more than by
glucose: **(E)-2-heptenal–GSH hemithioacetal M2 ≈ 2.85 mmol/L (4.8 % of GSH) at ~0.16 h**, decadienal
analogue N4 ≈ 0.9 mmol/L; these decay within 0.5 h into slower Cys-Gly thiazolidines (M1 ≈ 0.45 mM at
0.5 h, still 0.13 mM at 3 h). Hexanal is a weak competitor (L1 ≤ 0.03 mM). All intermediates are ≤ ~0.01
mM by 2 h. The quantification is GSH-equivalent with an assumed response, so the absolute scale is
semi-quantitative [!]; the timing and the ordering are the robust content.

---

## 2. METHODS [M] (pp. 2–3)

| item | as printed |
|---|---|
| Reactants | L-glutathione (99 %), D-glucose (98 %), hexanal (95 %), (E)-2-heptenal (95 %), (E,E)-2,4-decadienal (93 %), J&K |
| Target systems | Glu-GSH-Hex, Glu-GSH-Hept, Glu-GSH-Dec: **GSH 0.30 mmol + glucose 0.30 mmol + aldehyde 0.70 mmol** in **5 mL phosphate buffer (0.2 M, pH 6.5)** [Z: 60 / 60 / 140 mmol/L] |
| Controls | GSH-Glu, GSH-Hex, GSH-Hept, GSH-Dec (amounts "according to the respective systems") |
| Vessel | **15 mL pressure-resistant tube, sealed** (headspace ~10 mL air, not flushed) |
| Heating | **140 °C with stirring for 0.16, 0.25, 0.5, 1, 2, 3 or 4 h**; heating device and come-up not stated [!]; **two replicates** |
| pH | buffer label only; no measured pH [!]. Aldehydes at 140 mM are a separate phase (hexanal solubility ~5 g/L ≈ 50 mM; decadienal far lower) [Z, general knowledge] |
| Quantification (HPLC-MS SIM, [M+H]⁺) | GSH used as "internal standard". (1) 0.16 h mixtures diluted 0, 4, 10, 20, 30, 100-fold; slopes of peak area vs dilution give ki (intermediate) and k0 (GSH); ρ = ki/k0 "response factor". (2) GSH calibration factor f0 (mmol/L per area) from the 0.16 h mixture spiked with GSH. (3) **Ci = (ρ/n) × f0 × Ai**, n = number of cysteine residues in the intermediate |
| HPLC-MS | LCQ Deca XP MAX ion trap, ESI⁺; XBridge Amide 4.6 × 150 mm, ACN : 10 mM ammonium formate 75:25, 1 mL/min, 40 °C |
| Prepared standards | L1 (GSH 0.6 + hexanal 0.6 mmol, 90 °C 80 min), M2 (GSH 0.9 + heptenal 0.6 mmol, 90 °C 10 min), M3 (GSH 1.5 + heptenal 0.6 mmol, 90 °C 80 min), N4/N5 (GSH 1.5 + decadienal 1.5 mmol, 90 °C 5 min), each in 5 mL 0.2 M phosphate pH 6.5; purified on AG 50W-X4, HP 20SS or polyamide; purities (HPLC-ELSD) L1 99 %, M2 99 %, M3 94 %, N4 94 %, N5 94 %; yields 0.1, 0.1, 0.06, 0.02, 0.08 g |
| Volatiles from adducts | L1, M2, M3, N4, N5 at **0.06 M** in phosphate buffer, **140 °C, 2 h**; SPME DVB/CAR/PDMS 50 °C (10 + 20 min); DB-Wax; **relative peak area (% of total)** only |

**[!] The quantification is not a true calibration.** The ratio of two dilution slopes, ki/k0, equals
(response × concentration of the intermediate)/(response × concentration of GSH) in the same mixture; it
does not separate response from concentration. Combined with f0 and Ai it reduces to an assumption of
**equal molar ESI response per cysteine residue** (or the result is circular). Read every concentration
below as "GSH-equivalent, per cysteine", uncertain by an unknown factor (plausibly 2–10 × for sugar
adducts vs lipophilic adducts in ESI⁺) [Z]. Ratios **within one compound over time** are robust;
ratios **between compounds** carry the response uncertainty.

---

## 3. RESULTS

### 3.1 Identities [M] (Figs. 2–5, NMR in §2.3)

| code | system | structure (authors) | [M+H]⁺ |
|---|---|---|---|
| G1 | all with glucose | thiazolidine of **Cys-Gly** + C3 sugar fragment (glyceraldehyde-type) | 251 |
| G2 | all with glucose | thiazolidine of Cys-Gly + C4 sugar fragment | 281 |
| G3 | all with glucose | thiazolidine of Cys-Gly + intact glucose (reported before by Jerić & Horvat 2009) | 341 |
| G4 | all with glucose | **N-(1-deoxy-D-fructos-1-yl)-glutathione (GSH-Amadori)** | 470 |
| L1 | Hex | 2-pentylthiazolidine of Cys-Gly | 261 |
| M1 / N1 | Hept / Dec | 2-alkenylthiazolidine of Cys-Gly | 273 / 313 |
| N2 | Dec | = L1 structure (hexanal from decadienal retro-aldol) | 261 |
| M2 / N4 | Hept / Dec | **hemithioacetal of GSH-SH with the aldehyde carbonyl** (olefin retained, E) | 420 / 460 |
| M3 / N5 | Hept / Dec | complex thiazolidine: 1 GSH (Michael addition at C-3) + 1 Cys-Gly (thiazolidine) | 580 / 620 |
| N3 | Dec | 2 Cys-Gly | (MS only) |
| N6 | Dec | 1 GSH + 2 Cys-Gly | (MS only) |

GSSG is visible as a peak in the GSH-Glu and GSH-Dec control chromatograms (Fig. 1a, g) [Q]; not
quantified.

### 3.2 Fig. 1(h–j) (p. 4): concentrations, mmol/L, GSH-equivalent [F]

x-axis as plotted: the first non-zero point sits at **~0.10 h on the axis**, while Methods give 0.16 h as
the first time [!]; I label it 0.16 h. No 4 h point is plotted (axis ends at 3 h) [!]. Error bars only on
a few points. Read error: **(h) ±0.005 mM; (i) ±0.05 mM; (j) ±0.02 mM**.

**(h) Glu-GSH-Hex**

| t (h) | G1 | G2 | G3 | G4 (GSH-Amadori) | L1 |
|---|---|---|---|---|---|
| 0.16 | 0.003 | ~0.02 | 0.022 | **0.175** | 0.022 |
| 0.25 | 0.078 | 0.123 | **0.258** | 0.090 | 0.029 |
| 0.5 | **0.143** | **0.280** | 0.077 | 0.017 | 0.025 |
| 1 | 0.115 | 0.200 | 0.038 | 0.003 | 0.017 |
| 2 | ~0 | ~0.005 | ~0 | ~0 | 0.007 |
| 3 | ~0 | ~0 | ~0 | ~0 | 0.003 |

**(i) Glu-GSH-Hept**

| t (h) | M1 | M2 (hemithioacetal) | M3 | G1–G4 |
|---|---|---|---|---|
| 0.16 | ~0.02 | **2.85** | **0.60** | each ≤ 0.12 |
| 0.25 | 0.36 | 0.90 | 0.52 | G3 ~0.17, others ≤ 0.1 |
| 0.5 | **0.45** | ~0.10 | 0.10 | G2 ~0.25, others ≤ 0.1 |
| 1 | 0.35 | ~0.02 | ~0.02 | G2 ~0.12, others ≤ 0.05 |
| 2 | 0.22 | ~0 | ~0 | ~0 |
| 3 | 0.13 | ~0 | ~0 | ~0 |

**(j) Glu-GSH-Dec**

| t (h) | N4 (hemithioacetal) | N5 | N6 | N1, N2, N3 | G4 | G3 | G2 |
|---|---|---|---|---|---|---|---|
| 0.16 | **0.90** | 0.27 | 0.09 | ≤ 0.02 | 0.13 | ~0.02 | ~0.02 |
| 0.25 | 0.26 | ~0.15 | **0.45** | ≤ 0.06 | ~0.05 | 0.19 | ~0.06 |
| 0.5 | ~0.01 | ~0.05 | 0.41 | ≤ 0.09 (N3 ~0.085) | ~0.01 | ~0.08 | **0.20** |
| 1 | ~0 | ~0.02 | 0.12 | ≤ 0.03 (N2 peaks here per text) | ~0 | ~0.03 | ~0.09 |
| 2 | ~0 | ~0 | ~0.01 | ~0.01 | ~0 | ~0 | ~0.01 |

Peak-time order as printed in text [Q]: Hex system G4 → G3 → G1/G2; L1 at 0.25 h. Hept: M2/M3 at 0.16 h,
M1 at 0.5 h. Dec: N4/N5 at 0.16 h, N6 0.25 h, N1/N3, N2 at 1 h.

**Derived [Z]** (charge 60 mmol/L GSH):

| quantity | value |
|---|---|
| GSH-Amadori peak / GSH charged (Hex system) | 0.175/60 = **0.29 %** |
| Σ Cys-Gly–sugar thiazolidines G1+G2+G3 at 0.5 h (Hex) | 0.50 mM = 0.8 % |
| M2 + M3 at 0.16 h (Hept), GSH-equivalent | 3.45 mM = **5.8 %** of GSH |
| N4 + N5 + N6 at 0.16 h (Dec) | 1.26 mM = 2.1 % |
| G4 apparent decay, Hex, 0.16 → 0.5 h | k ≈ ln(0.175/0.017)/0.34 h ≈ 6.9 h⁻¹ ≈ **0.11 min⁻¹** |
| M2 apparent decay, 0.16 → 0.25 h | ≈ ln(2.85/0.90)/0.09 h ≈ 13 h⁻¹ ≈ 0.2 min⁻¹ (two points, ±50 %) |
| M1 decay, 0.5 → 3 h | ln(0.45/0.13)/2.5 h ≈ 0.50 h⁻¹ ≈ 8 × 10⁻³ min⁻¹ |

The G4 decay rate (~0.1 min⁻¹ at 140 °C) agrees with the GSH-Amadori decay in du2025 at 145 °C (falls
~3.5 × per 10 min, ≈ 0.12–0.18 min⁻¹) [Z].

### 3.3 Volatiles from re-heated adducts (§3.4) [Q]

Table S1 is not on disk. Text: heated L1, M2, M3, N4, N5 (0.06 M, 140 °C, 2 h) released mainly **the
original aldehydes** (hexanal, (E)-2-heptenal, (E,E)-2,4-decadienal; (E)-2-heptenoic acid), i.e. the
adducts are **reversible**; plus retro-aldol/oxidation products and alkyl S/N compounds (2-hexylthiophene,
2-ethylpyridine, 1-hexanethiol, 1-(2-thienyl)-1-propanone/-hexanone, six 2,5-dialkylthiophenes such as
2-ethyl-5-(1-methylbutyl)thiophene). The same alkyl compounds appear in GSH + aldehyde without glucose
("data not shown"). No MFT or FFT is reported from the adducts.

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Model context.** `GSH` → H₂S + fragments only (`r_gsh_h2s`, Zheng & Ho 1994, 80–110 °C); no GSH
glycation, no γ-glutamyl cleavage, no thiol–aldehyde adducts; the lipid lane's aldehydes do not react with
any thiol. Nothing here was run against the model.

**(1) Glutathione's fate in a sugar pot (qualitative-to-semi-quantitative, usable as a contradiction).**
At 140 °C the GSH-Amadori compound is a **minor (≤ 0.3 %), short-lived (τ ≈ 10 min)** pool; the
Cys-Gly products (thiazolidines) and Cys-Gly itself appear within 10 min. Glutathione's sugar chemistry
is therefore dominated by **γ-glutamyl loss → Cys-Gly → thiazolidine (TTCA-like) chemistry**, not by a
stable GSH-Amadori. This parallels the free-cysteine TTCA/Cys-Amadori chemistry the model does carry
(`cao2017_extraction.md`, `zhang2025b_extraction.md`). A minimal GSH extension would be
GSH → Cys-Gly (+ pGlu), with Cys-Gly then treated like cysteine; the rate is constrained by
`zhou2023c_extraction.md` (§4), not by this paper (no free GSH or Cys-Gly time course here).

**(2) Lipid × sulfur: thiol capture by α,β-unsaturated aldehydes (semi-quantitative).** 140 mM
(E)-2-heptenal binds ~6 % of 60 mM GSH within ~10 min as hemithioacetal/Michael adducts — **≈ 20 × the
GSH-Amadori pool**; hexanal (saturated) barely competes. The adducts are reversible (re-heating releases
the aldehyde). For a model coupling of lipid aldehydes to thiols, this supplies: (i) the **ordering**
2-alkenal ≈ 2,4-dienal ≫ alkanal for thiol capture; (ii) fast pre-equilibrium (≤ 10 min at 140 °C);
(iii) an order of magnitude for the bound fraction at a 2.3 : 1 aldehyde : GSH excess. It cannot give an
equilibrium constant (aldehyde phase behaviour unknown, response factors assumed).

**(3) Not usable.** Absolute concentrations as fit targets (response assumption, two replicates, a
plotted first time that disagrees with Methods); any volatile amount (relative %; Table S1 absent); the
pathway schemes (Fig. 6) as rates; nothing on H₂S, MFT or FFT.

**Caveats, collected.**
1. Buffer label pH 6.5; 60 mM GSH + 140 mM aldehyde, actual pH unknown.
2. Aldehydes far above aqueous solubility at 140 mM; two-phase pots.
3. Quantification = GSH-equivalent per cysteine; between-compound ratios uncertain.
4. Free GSH never quantified over time (used as calibrant); GSSG seen but not quantified.
5. First plotted time ~0.10 h vs 0.16 h in Methods; 4 h point not plotted.
