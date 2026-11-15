# Morales & Jiménez-Pérez 2001 — EXTRACTION
### *"Free radical scavenging capacity of Maillard reaction products as related to colour and fluorescence"*
### Glucose or lactose (0.1 M) alone or with alanine, glycine or lysine (0.1 M) in 100 mM phosphate, started at pH 7 and **not controlled**, stoppered Pyrex flasks in a **100 °C** glycerol bath for **3, 6, 12, 24 h**; pH, A420, CIELab, fluorescence (ex 347 / em 415 nm, quinine-sulphate mg/L) and **DPPH· scavenging** (200 µL sample + 1 mL DPPH· 74 mg/L in chloroform, 1 h, A520) → TC50, the heating time at which the heated mixture removes 50 % of the DPPH·. **All results are figures; no table.** Lysine MRPs reach TC50 in ~3.5 h, alanine/glycine ~4.7–6.7 h, sugars alone ~11.5 h; glucose vs lactose differs by ≤ 0.15 dex in TC50. **No lipid in any pot.** Read 2026-11-14.

**Full read of `data/articles/morales2001.pdf`** (7 pp., Food Chem. 72 (2001) 119–125; Elsevier PDF with
a text layer whose ligatures are broken: "¯" = fl, "®" = fi, "" = ff, "±" = en dash, "ÿ" = minus).
Prose from `pdftotext -layout`. Methods page re-read from a 200 dpi render (p. 120): the sample
volumes the text layer prints as "100 ml" and "200 ml" are **100 µl** and **200 µl**. Every number in
§3 was read from 170 dpi renders of pp. 121, 123 and 124 (Figs. 1, 2, 4, 5, 6); Fig. 3 (CIELab) is
described only qualitatively.

**Provenance:** **[M]** measured, printed · **[Q]** qualitative · **[Z]** derived by me (incl. values
read off figures; reading uncertainty stated per figure) · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Francisco J. Morales,* Salvio Jiménez-Pérez |
| Affiliation | Instituto del Frío (C.S.I.C.), Ciudad Universitaria, Madrid |
| Title | "Free radical scavenging capacity of Maillard reaction products as related to colour and fluorescence" |
| Venue | *Food Chemistry* **72 (2001) 119–125** (Analytical, Nutritional and Clinical Methods Section) |
| DOI | **10.1016/S0308-8146(00)00239-9** (from the printed PII S0308-8146(00)00239-9; matches the expected DOI) |
| Dates | received 21 March 2000; revised 15 June 2000; accepted 15 June 2000 |
| Funding | EU FAIR-CT-96-1080 ("Optimisation of the Maillard Reaction") |
| Registry | not in `data/keys/papers.yml` (2026-11-14) |

---

## 1. THE ONE-PARAGRAPH ANSWER

Equimolar sugar/amino-acid solutions (0.1 M each, 100 mM phosphate, initial pH 7, 100 °C, 3–24 h, pH
allowed to drift) develop DPPH·-scavenging activity fast: in every sugar + amino-acid pot the DPPH·
remaining falls from ~94 % (unheated) to ≤ 17 % by 12 h and plateaus at ~8–12 % by 24 h. Ranked by
TC50 (h of heating to reach 50 % DPPH· remaining): **glucose-lysine 3.4 ≈ lactose-lysine 3.7 <
glucose-glycine 4.7 ≈ glucose-alanine 4.7 < lactose-alanine 5.6 < lactose-glycine 6.7 ≪ glucose 11.4 ≈
lactose 11.7**. Lysine pots also acidify most (pH 4.33 GL, 4.11 LL at 24 h) and brown most. Scavenging
tracks fluorescence (which plateaus after ~6 h), not browning (which keeps rising roughly linearly).
The authors' point is analytical: fluorescence is a better proxy than A420 for antiradical MRPs. The
paper measures an antiradical assay in chloroform, **not lipid oxidation**, and varies sugar (glucose
vs lactose) and amino acid (Ala/Gly/Lys), **not ribose, cysteine, thiamine or a controlled pH**.

---

## 2. METHODS [M] (pp. 120–121)

| item | as printed |
|---|---|
| Chemicals | D-(+)-glucose, L-lysine monohydrochloride, L-alanine, L-glycine "monohydrochloride", DPPH· (Sigma); D-(+)-lactose monohydrate, chloroform (Panreac) |
| Pots | **sugar 0.1 mol + amino acid 0.1 mol in 1 L of 100 mM phosphate buffer (pH 7)** → 0.1 M each, 1:1. Eight systems: G, GA, GG, GL, L, LA, LG, LL (sugar alone = G, L) |
| pH control | **none**: "heated without pH control"; pH followed over time (Fig. 1) |
| Heating | tightly stoppered Pyrex flasks immersed in a **glycerol bath at 100 °C**; sampled at **3, 6, 12, 24 h**, cooled in ice. Come-up not stated; headspace/atmosphere not stated (air) |
| Browning | A420, Shimadzu UV-1601, 1 cm cell, "after appropriate dilution" (dilution not stated) |
| Colour | MiniScan MS/S-4000S tristimulus, CIELab, D65 / 10°, 5 cm Petri dish; E = (L*² + a*² + b*²)^½, C = (a*² + b*²)^½; mean of four determinations at 22–24 °C |
| Fluorescence | **100 µL** sample diluted with 10 mL Milli-Q water; ex 347 / em 415 nm; Kontron SMF-25; 1 cm quartz; calibrated vs quinine sulphate in 0.1 M H2SO4, results as **quinine sulphate mg/L**; mean of three readings |
| DPPH· assay | after Cämmerer & Kroh (1996), modified: **200 µL** sample + **1 mL DPPH· 74 mg/L in chloroform** (initial A520 ≈ 1.8); shaken vigorously **1 h** (steady state; LL-24 h reached 50 % in ~8.5 min, plateau 20 min; all others < 60 min); A520 of the organic layer at **30 °C**; [DPPH·]t (mg/100 mL) = 0.2563·A520 + 0.0203 (r² = 0.9995); %DPPH·rem = [DPPH·]t/[DPPH·]0 × 100 |
| TC50 | heating time (h) needed for the heated mixture to bring DPPH· to 50 % remaining, read from %DPPH·rem vs heating time |
| Interference checks | coloured MRPs did not colour the chloroform layer after 60 min |
| Repeatability | highest CV 3.2 ± 0.05 % (n = 6) for LG 6 h (DPPH·) |
| n / stats | "at least in duplicate"; ANOVA + Duncan, 95 % |

---

## 3. RESULTS [M] (text) and [Z] (read from figures)

### 3.1 Printed numbers [M]
- pH after 24 h: **LL 4.11, GL 4.33** (the lowest; lysine systems significantly lower than the others).
- DPPH· kinetics: LL-24 h sample reaches 50 % DPPH· remaining at ~8.5 min of assay, plateau 20 min.
- CV of the DPPH· assay ≤ 3.2 % (n = 6, LG 6 h).
- Antiradical groups (from Fig. 6, as stated): GL, LL high; GA, GG medium; LA, LG medium–low; G, L low.
- Chroma maximum at 3 h for LL, GL, GG; at 6 h for GA, LG, LA.

### 3.2 Fig. 1 — pH vs heating time [Z] (read ±0.05 pH; symbols: empty = glucose, full = lactose)

| system | 0 h | 3 h | 6 h | 12 h | 24 h |
|---|---:|---:|---:|---:|---:|
| G | 7.05 | 6.97 | 6.83 | 6.67 | 6.21 |
| L | 7.00 | 6.90 | 6.68 | 6.33 | 5.75 |
| GA | 7.00 | 6.75 | 6.42 | 6.12 | 5.45 |
| LA | 7.00 | 6.80 | 6.52 | 6.18 | 5.56 |
| GG | 7.00 | 6.66 | 6.25 | 5.84 | 5.25 |
| LG | 7.00 | 6.70 | 6.33 | 5.84 | 5.23 |
| GL | 6.81 | 6.08 | 5.37 | 4.75 | 4.34 (printed 4.33) |
| LL | 6.77 | 6.08 | 5.33 | 4.57 | 4.11 (printed 4.11) |

Lysine pots start at ~6.8, not 7.0 (lysine·HCl in 100 mM phosphate) [!].

### 3.3 Fig. 2 — browning, A420 [Z] (read ±0.005; dilution factor not stated)

| system | 3 h | 6 h | 12 h | 24 h |
|---|---:|---:|---:|---:|
| G | 0.000 | 0.005 | 0.018 | 0.045 |
| L | 0.000 | 0.008 | 0.020 | 0.053 |
| GA | 0.020 | 0.055 | 0.115 | 0.200 |
| LA | 0.005 | 0.030 | 0.068 | 0.160 |
| GG | 0.030 | 0.085 | 0.153 | 0.256 |
| LG | 0.015 | 0.065 | 0.125 | 0.235 |
| GL | 0.100 | 0.220 | 0.340 | 0.487 |
| LL | 0.045 | 0.115 | 0.200 | 0.335 |

### 3.4 Fig. 4 — fluorescence, quinine-sulphate mg/L [Z] (read ±5 mg/L)

| system | 3 h | 6 h | 12 h | 24 h |
|---|---:|---:|---:|---:|
| G | 2 | 3 | 5 | 12 |
| L | 4 | 10 | 20 | 61 |
| GA | 55 | 109 | 124 | 130 |
| LA | 26 | 82 | 122 | 135 |
| GG | 77 | 141 | 172 | 179 |
| LG | 72 | 155 | 189 | 198 |
| GL | 178 | 215 | 214 | 195 |
| LL | 155 | 241 | 245 | 223 |

### 3.5 Fig. 5 — % DPPH· remaining after 1 h of assay, vs heating time [Z] (read ±2 %; (A) glucose systems, (B) lactose systems; sugar-alone curves are not plotted)

| system | 0 h | 3 h | 6 h | 12 h | 24 h |
|---|---:|---:|---:|---:|---:|
| GA | 94 | 65 | 43.5 | 17 | 11.5 |
| GG | 94 | 67.5 | 42 | 13 | 8 |
| GL | 94 | 56 | 29 | 11.5 | 9 |
| LA | 94 | 85.5 | 48 | 17 | 7.5 |
| LG | 93.5 | 88.5 | 57.5 | 14 | 9.5 |
| LL | 95.5 | 59 | 25 | 9.5 | 8 |

### 3.6 Fig. 6 — TC50 (h), bar ± error bar [Z] (read ±0.1 h; the error bar's meaning — SD, SE, n — is not stated [!])

| system | TC50 (h) | bar |
|---|---:|---:|
| G | 11.4 | ±0.3 |
| GA | 4.7 | ±0.2 |
| GG | 4.7 | ±0.15 |
| GL | 3.4 | ±0.15 |
| L | 11.7 | ±0.3 |
| LA | 5.6 | ±0.2 |
| LG | 6.7 | ±0.3 |
| LL | 3.7 | ±0.2 |

### 3.7 Fig. 3 — CIELab [Q]
E falls (darkening) in all pots; GL, GG, LL darken faster than LG, GA, LA but converge by 24 h. a*b*
trajectories loop from yellow-brown to orange-brown and, for GL and LL at long times, purplish red.

---

## 4. DERIVED [Z]

### 4.1 log10 contrasts in TC50 (positive = the first-named pot is a faster / stronger scavenger former)

| contrast | log10 ratio |
|---|---:|
| amino acid present vs sugar alone: G vs GL | +0.53 |
| G vs GA / GG | +0.38 / +0.38 |
| L vs LL | +0.50 |
| L vs LA / LG | +0.32 / +0.24 |
| **sugar identity (glucose vs lactose), fixed amino acid: Lys** | **+0.04** |
| Ala | +0.08 |
| Gly | +0.15 |
| sugar alone G vs L | +0.01 |
| amino acid identity, glucose: Lys vs Ala/Gly | +0.14 |
| amino acid identity, lactose: Lys vs Ala / Gly | +0.18 / +0.26 |

At 6 h of heating, % DPPH· remaining glucose vs lactose: Lys 29 vs 25 (−0.06 dex), Ala 43.5 vs 48
(+0.04), Gly 42 vs 57.5 (+0.14). By 12–24 h every sugar + amino-acid pot sits at 8–17 %: the
end-points converge regardless of sugar or amino acid.

### 4.2 Thermal-dose comparison with the P3 pot (assumption-labelled)
With an **assumed** apparent activation energy Ea for "antiradical MRP formation" of 80 / 100 / 120
kJ/mol (not measured in this paper or, to my knowledge, anywhere in the corpus), 20 min at 140 °C equals
**~4 / ~8 / ~14 h at 100 °C** (factor exp[Ea/R·(1/373.15 − 1/413.15)] = 12 / 23 / 42). The P3 cook
therefore probably lies **at or past TC50** of these pots — in the regime where sugar identity
(glucose vs lactose) barely matters and the scavenging endpoint has converged. This is a declared
assumption, not a measurement.

---

## 5. USE IN THIS REPO

**Question it was downloaded for (L1 scoping):** do Maillard reaction products act as antioxidants
on lipid oxidation in a heated pot, and which sugar / amino acid / T / time / pH controls that?

- **What it measures is not lipid oxidation.** DPPH· in chloroform is a radical-quenching capacity
  assay of the aqueous MRP mixture. No lipid, no hydroperoxide, no aldehyde, no headspace. It cannot
  give a ratio of hexanal, nonanal, 2-pentylfuran or decadienal with vs without MRPs.
- **Sugar lever (P3: glucose vs ribose):** ribose is not tested. The only sugar contrast is glucose
  vs lactose, worth **≤ 0.15 dex in TC50** and converging to the same plateau. A mechanism pointer
  that antiradical capacity is driven by the amine partner (lysine ≫ Ala/Gly, +0.14 to +0.26 dex) and
  by extent of reaction, more than by sugar identity. A pentose would react faster than glucose, which
  would only move the pot further onto the plateau [Q, not measured here].
- **Sulfur lever (cysteine vs thiamine):** not tested. Nothing.
- **pH lever:** not tested as a variable; pH drifts from 7 to 4.1–6.2. Nothing usable on pH 6 vs 7.
- **Distance from the P3 pot:** 100 °C / 3–24 h vs 140 °C / 20 min (thermally overlapping only under
  the assumed Ea of §4.2); free amino acids at 0.1 M vs a 15 % pea-isolate protein matrix whose lysine
  residues (~0.07 M [Z], assuming ~7 g Lys/100 g protein) dominate the amine pool regardless of the
  additive lever; no lipid; no sulfur; unbuffered drift vs held pH.
- **Usable for a declared band?** **No.** It supports only the qualitative statement in the L1
  scoping note's direction: MRP antiradical capacity forms within hours at 100 °C from any reducing
  sugar + amine, so in a heated protein + sugar pot an MRP antioxidant term is plausible and is
  insensitive to glucose-vs-lactose sugar identity (≤ 0.15 dex on its own assay scale). Translating a
  DPPH· TC50 shift into a hexanal ratio would need a lipid pot that this paper does not have.
- **Hold-out discipline [!].** `farmer1992` and `elmore2002` are designated hold-outs; do not combine
  them with this paper into any proposed band. (No overlap of lab: Madrid, CSIC.)

---

## 6. FLAGS

1. **Text-layer ligature damage** throughout ("¯uorescence", "eect"); µ rendered as "m" in the sample
   volumes (100 µl / 200 µl, verified on the image).
2. **"glucose-glycine (LG)"** in Methods 2.2 is a typo for lactose-glycine (all figures and the
   results treat LG as lactose-glycine).
3. **"significant (P > 0.05)"** (twice, p. 121) — inverted inequality; read as P < 0.05.
4. **"previously heated at different temperatures"** (2.3.4) — only 100 °C is reported anywhere.
5. **Figure-only data.** Every number in §3.2–3.6 is my reading of a scanned figure; uncertainty is
   the stated reading error, not experimental error. Fig. 6 error bars are undefined.
6. **Sugar-alone DPPH· curves are not plotted** (Fig. 5 shows six curves) though G and L TC50 values
   appear in Fig. 6.
7. **0 h DPPH· remaining ≈ 94 %**, not 100 %: the unheated mixtures (or the aqueous aliquot itself)
   remove ~5 % of DPPH·.
8. **Lysine·HCl** lowers the starting pH (~6.8); part of the lysine pots' low TC50 may come from the
   different pH path, not only the extra ε-amine.
