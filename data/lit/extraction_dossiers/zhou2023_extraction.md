# Zhou, Xia, Cui, Zhai, Zhang, Hayat, Zhang & Ho 2023 (10.1021/acs.jafc.2c08360) — Wave K3 extraction 2026-08-28

**THE CORRECT PAPER HAS ARRIVED.** Source PDF: `data/articles/Zhou2023b.pdf` (10 pp.)
**plus** `data/articles/zhou2023b_supporting.pdf` (9 pp., the Supporting Information).

> Zhou, T.; Xia, X.; Cui, H.; Zhai, Y.; Zhang, F.; Hayat, K.; Zhang, X.\*; Ho, C.-T.\*
> "**Cysteine-Induced pH-Dependent Formation of Thiols and Sulfides or 2-Acetylthiazole and
> Pyrazines during Thermal Treatment of N-(1-Deoxy-D-xylulos-1-yl)-alanine**"
> *J. Agric. Food Chem.* **2023**, 71 (5), 2472−2481. Jiangnan University (Wuxi) + Rutgers.
> Received Nov 28 2022; revised Jan 10 2023; accepted Jan 12 2023; published Jan 25 2023.
> Downloaded stamp on the PDF: "by EUROPEAN COMMISSION JRC user on 28 August 2026".

**Read method.** Full text layer read directly (`k3_zhou2023b.txt`); Table 1 and Tables S1/S2
came out of the text layer intact and are transcribed verbatim. Figures 2 and 3 have **no
numeric text layer** and were digitised from **300 dpi page renders** (`k3img/zhou_p-05.png`,
`k3img/zhou_p-06.png`, crops `z_fig2_zoom.png`, `z_fig3_full.png`, `z_fig3_bot.png`).
Every figure-derived number below is tagged **[fig]**.

> ⚠️ **FILENAME NOTE.** The previous `Zhou2023_extraction.md` in this scratchpad documented a
> *different* paper (`data/articles/Zhou2023.pdf` = **Zhu et al., kafirin × pickle-like odorants,
> Food Chem. 400 (2023) 133854**), which was mis-filed under this DOI. That dossier has been
> **preserved unchanged** as `zhu2023_kafirin_binding_extraction.md` — it holds the only
> MFT–protein *binding constants* in the corpus (K ≈ 1.1 × 10³ M⁻¹ at 25 °C) and must not be
> lost. `data/articles/Zhou2023.pdf` is still the wrong file for DOI 10.1021/acs.jafc.2c08360.

---

## 0. VERDICT UP FRONT

**This is the paper the sulfur module was missing, and it delivers on all three counts the
brief asked for.**

1. **A purified, single-compound Amadori feedstock.** Not "sugar + amino acid heated together":
   a chromatographically purified, NMR- and MS/MS-confirmed **N-(1-deoxy-D-xylulos-1-yl)-alanine**
   at a known 20 mmol/L, heated with or without 20 mmol/L Cys. This is a **step-level fed-intermediate
   system**: the model can be initialised *at the Amadori node* with zero free sugar and zero
   ambiguity about how much ARP was present.
2. **A complete pH 6 / 7 / 8 series for 22 quantified volatiles**, both with and without Cys,
   at one fixed 120 °C / 60 min (Table 1, µg/L). Six condition columns.
3. **A second, deeper fed-intermediate system**: purified **Cys 20 mmol/L + MGO 20 mmol/L, 1:1**,
   pH 6/7/8, **0 / 30 / 60 / 90 min at 120 °C** — a genuine *concentration–time* series at three
   pH values for 8 compounds (Figure 3). This is the only α-dicarbonyl-fed, time-resolved,
   pH-resolved dataset in the whole corpus.

**Plus, from the SI: an odor-threshold list in water (µg/L) with 13 entries** (Table S2),
including MFT = **0.005 µg/L** and bis(2-methyl-3-furyl) disulfide = **0.00032 µg/L**.

**★ THE HEADLINE RESULT — MFT IS NON-MONOTONE IN pH WITH A MAXIMUM AT pH 7**, and this is a
**direction conflict with Hofmann & Schieberle 1998** over the pH region where the two overlap.
See §6. Do not paper over it: it is a real, informative structural constraint.

**No rate constants. No activation energies. No Arrhenius data. One temperature (120 °C)
throughout.** Figure 3 is a time series but the paper fits nothing to it; any rate extracted
from it is **your** fit, not theirs, and must be labelled as such.

---

## 1. SYSTEM DEFINITIONS — verbatim

### 1.1 Preparation and purification of the ARP (Methods, p. 2473)

> "Xyl (**6.0 g**) and Ala (**1.78 g**) were dissolved in **100 mL of deionized water** with a
> **molar ratio of 2:1**, and the pH of the mixture was adjusted to **7.5** with NaOH (3 mol/L).
> After the above operation, the mixture was heated under normal pressure for **60 min at 80 °C**.
> Then the reaction liquid was dehydrated with a rotary evaporator under vacuum. After
> dehydration at **80 °C for 20 min**, the reaction was terminated immediately by placing it in an
> ice water bath and the obtained mixture was redissolved in **50 mL** of ultrapure water for
> further purification."

> "Before purification, the pH of the system containing the ARP was adjusted to **4** with
> 1 mol/L HCl and then loaded onto a column filled with **H⁺-form Dowex 50WX8** ion-exchange resin.
> After removing sugar with deionized water (2 mL/min), ARP was eluted from the mixture using
> **0.1 mol/L ammonium hydroxide** at a flow rate of 1 mL/min. Finally, the solution containing
> only purified ARP was gathered through HPLC detection and **lyophilized into powder**. The
> purified ARP was identified as **N-(1-deoxy-D-xylulos-1-yl)-alanine** through UPLC-Q-TOF/MS
> and NMR."

Confirmation data (SI): MS/MS parent ion **[M+H]⁺ m/z 222.0911**; fragments 204.0890 (−H₂O),
186.0781 (−2H₂O), 168.0673 (−3H₂O), 176.0968 (−H₂O−CO), 158.0876, 140.0758, 131.0627
(−H₂O−COOH), 102.0590 and 90.0584 (alanine residue). ¹H/¹³C NMR shifts for **keto-form,
α-anomer and β-anomer** are all given (SI Fig. S2c is the isomerisation scheme) — i.e. the
feedstock is a **three-species tautomeric equilibrium**, not a single structure. Keto-form
¹³C C-2 at **208.75 ppm**; α-anomer C-2 at **103.37 ppm**; β-anomer C-2 at **105.45 ppm**.

### 1.2 SYSTEM A — heat treatment of ARP ± Cys (Methods, p. 2473) — the Table 1 system

> "The aqueous solutions of **ARP (20 mmol/L)** with the presence or the absence of added
> **Cys (20 mmol/L)** were prepared, and the pH levels of the systems were adjusted to different
> pH values (**6, 7, and 8**) using NaOH (3 mol/L). The reaction was carried out at **120 °C for
> 1 h**, and then the reacted mixture was cooled in an ice bath..."

- Solvent: **deionized water**. ⚠️ **NO BUFFER.** pH is *initial* only and drifts hard (see §3).
- Molar ratio ARP:Cys = **1:1**.
- Six conditions: {pH 6, 7, 8} × {ARP alone, ARP + Cys}.
- n = 3 for GC-MS ("Each experiment was repeated three times... mean values ± standard
  deviation... p < 0.05", Statistical Analysis, p. 2474). n = 4 parallel samples for the
  electronic nose.

### 1.3 SYSTEM B — Cys + MGO model (Methods, p. 2473) — the Figure 3 system

> "**Cys (20 mmol/L)** and **MGO (20 mmol/L)** were dissolved in deionized water with a
> **molar ratio of 1:1**, and the pH of the solution was adjusted to **6, 7, and 8** using NaOH
> (3 mol/L). The reactions of the three mixtures were all performed for **30, 60, and 90 min at
> 120 °C**, respectively, and then immediately cooled in an ice bath after completion of the
> reaction."

- MGO supplied as **40 % aqueous solution** (J&K Scientific). Again **unbuffered**.
- Time points plotted: **0, 30, 60, 90 min** (t = 0 is drawn at zero concentration for every
  compound in every panel — it is an assumed origin, not a measured point).

### 1.4 Quantification method — read this before converting anything

> "HS−SPME was performed with a **carboxen/poly(dimethylsiloxane)/divinylbenzene
> (CAR/PDMS/DVB) fiber**... The reacted solution sample (**3 mL**) and saturated sodium chloride
> solution (**2 mL**) were mixed in a headspace vial. **1,2-Dichlorobenzene (0.0018 µg/µL in
> methanol) was added as an internal standard.** ... headspace extraction was performed at
> **60 °C for 30 min** under stirring at **300 rpm**. After sampling, the fiber was attached to
> the GC injector at 250 °C for 10 min, and the GC was run in the splitless mode."

Column: **DB-WAX 30 m × 0.25 mm × 0.25 µm**. He at 1.8 mL/min. Oven: 40 °C (3 min) → 5 °C/min
to 80 °C → 10 °C/min to 160 °C (0.5 min) → 2 °C/min to 175 °C → 10 °C/min to 230 °C (7 min).
EI 70 eV, m/z 35–450. Identification: **MS + RI + authentic standard (S)** for all 22 compounds
(Table 1 col. "identification methods" = "MS, RI, S" on every row).

Basis:
> "Different concentrations of the chemical standards and the fixed concentration of internal
> standards were mixed to construct standard curves. The **x-axis and the y-axis represent the
> concentration of the flavor compounds and the peak area ratio of the corresponding compounds
> to the internal standard**, respectively."

⚠️ **BASIS = µg per litre of the reacted LIQUID**, back-calculated through an external calibration
curve run under identical HS-SPME conditions (Table S1). It is **NOT** SIDA. Matrix-matched only
insofar as the standards were run in the same vial geometry. **Not method-commensurable with
Hofmann & Schieberle's SIDA numbers** — comparisons across the two must be made on *ratios and
directions*, not on absolute ppb.

---

## 2. ★ TABLE 1 (pp. 2474–2475) — THE FULL pH 6/7/8 SERIES, ALL 22 COMPOUNDS

**"Volatile Flavor Compounds of Thermally Treated ARP with or without Cys at 120 °C for 1 h under
Different Initial pH Values (6, 7, and 8)". Concentration in µg/L. mean ± SD, n = 3.**
Superscript letters = Duncan's multiple comparison across the row, p < 0.05. ND = not detected.
RI = linear retention index on DB-WAX (C7–C30 n-alkanes); KI = NIST/literature value.

| compound | RI | KI | pH 6 ARP | pH 6 ARP+Cys | pH 7 ARP | pH 7 ARP+Cys | pH 8 ARP | pH 8 ARP+Cys |
|---|---|---|---|---|---|---|---|---|
| 2,3-butanedione | 979 | 995 | 8.93 ± 0.53 ᵇ | **ND** | 33.63 ± 1.37 ᶜ | **ND** | 126.18 ± 4.27 ᶜ | 20.00 ± 1.55 ᵈ |
| **Furans** | | | | | | | | |
| furan | 797 | 802 | ND | ND | ND | 30.77 ± 0.97 ᵇ | 34.88 ± 1.75 ᵇ | 63.33 ± 3.44 ᶜ |
| 2-furfural | 1461 | 1460 | **1606.92 ± 41.62** ᵉ | 608.41 ± 22.97 ᵈ | 1339.37 ± 83.04 ᶜ | 473.07 ± 14.04 ᶜ | 436.63 ± 19.93 ᶜ | 224.98 ± 23.14 ᵇ |
| 2-acetylfuran | 1499 | 1501 | 32.34 ± 1.49 ᶜ | 6.67 ± 0.63 ᵇ | 62.17 ± 2.14 ᵈ | 6.63 ± 0.65 ᵇ | 93.48 ± 2.91 ᵉ | 9.81 ± 0.34 ᶜ |
| **Pyrazines** | | | | | | | | |
| pyrazine | 1205 | 1210 | ND | ND | ND | 6.06 ± 0.81 ᵇ | 3.99 ± 0.76 ᵇ | **177.05 ± 8.08** ᶜ |
| methylpyrazine | 1265 | 1263 | ND | ND | ND | 1.03 ± 0.08 ᵇ | 3.97 ± 0.15 ᶜ | **99.00 ± 1.96** ᶜ |
| 2,5-dimethylpyrazine | 1321 | 1315 | ND | ND | ND | ND | 0.272 ± 0.06 | 21.44 ± 2.00 |
| 2,6-dimethylpyrazine | 1328 | 1321 | ND | ND | ND | ND | ND | 0.794 ± 0.09 |
| **Thiophenes** | | | | | | | | |
| thiophene | 1014 | 1038 | ND | 13.73 ± 0.56 ᵇ | ND | 21.77 ± 1.64 ᶜ | ND | 16.95 ± 3.82 ᵇ |
| 2-methylthiophene | 1088 | 1117 | ND | 13.44 ± 0.10 ᵇ | ND | 56.28 ± 3.74 ᶜ | ND | 116.24 ± 10.91 ᶜ |
| 2-thiophene-carboxaldehyde | 1693 | 1677 | ND | 21.66 ± 1.43 ᵇ | ND | 45.59 ± 2.07 ᶜ | ND | 43.63 ± 2.62 ᶜ |
| 5-methyl-2-thiophene-carboxaldehyde | 1711 | 1759 | ND | 1.86 ± 0.35 ᵇ | ND | 2.92 ± 0.16 ᵇ | ND | 55.44 ± 8.39 ᶜ |
| thieno[3,2-b]thiophene | 1879 | 1843 | ND | 4.78 ± 0.83 ᶜ | ND | 5.47 ± 1.11 ᶜ | ND | 2.03 ± 0.66 ᵇ |
| **Thiazoles** | | | | | | | | |
| 2-methylthiazole | 1237 | 1239 | ND | 0.03 ± 0.01 ᵇ | ND | 0.21 ± 0.02 ᶜ | ND | 3.23 ± 0.14 ᶜ |
| thiazole | 1247 | 1259 | ND | 5.27 ± 0.22 ᵇ | ND | 10.90 ± 0.99 ᵇ | ND | **107.38 ± 8.50** ᶜ |
| 2,4-dimethylthiazole | 1315 | 1283 | ND | ND | ND | ND | ND | 0.77 ± 0.09 |
| **2-acetylthiazole** | 1647 | 1652 | ND | 2.43 ± 0.08 ᵇ | ND | 11.70 ± 2.14 ᵇ | ND | **582.34 ± 14.57** ᶜ |
| **Thiols and sulfides** | | | | | | | | |
| **2-methyl-3-furanthiol (MFT)** | 1303 | 1316 | ND | **696.99 ± 17.44** ᶜ | ND | **1588.57 ± 21.24** ᶜ | ND | **525.62 ± 16.88** ᵇ |
| **2-furfurylthiol (FFT)** | 1430 | 1430 | ND | **813.65 ± 13.37** ᶜ | ND | **757.965 ± 13.03** ᶜ | ND | **325.22 ± 14.61** ᵇ |
| di-2-furfuryl sulfide | 2212 | 2223 | ND | 6.55 ± 0.69 | ND | 8.68 ± 0.91 | ND | **ND** |
| bis(2-methyl-3-furyl) disulfide | 2148 | 2137 | ND | 59.70 ± 2.75 ᶜ | ND | 102.59 ± 4.17 ᶜ | ND | 50.07 ± 2.15 ᵇ |
| bis(2-furfuryl) disulfide | 2465 | 2458 | ND | 35.31 ± 1.63 | ND | 17.68 ± 0.76 | ND | **ND** |

Footnote a (verbatim): "Results are displayed as means ± standard deviation, data within a row
with different superscript letters denoted significantly different (p < 0.05) using Duncan's
multiple comparison test (n = 3)."

### 2.1 Derived quantities (MY arithmetic from Table 1 — clearly marked FITTED/DERIVED)

**Molar yields on the 20 mmol/L ARP basis** (mol product / mol ARP × 100). MWs used:
MFT & FFT 114.17 (both C₅H₆OS); 2-furfural 96.08; 2-acetylthiazole 127.16;
bis(2-methyl-3-furyl) disulfide **226.31** and bis(2-furfuryl) disulfide **226.31** (both
C₁₀H₁₀O₂S₂ — the two disulfides are isomers and share a molar mass);
di-2-furfuryl sulfide **194.25** (C₁₀H₁₀O₂S).

| quantity | pH 6 | pH 7 | pH 8 |
|---|---:|---:|---:|
| MFT, µmol/L | 6.105 | 13.914 | 4.604 |
| MFT, **mol% of ARP** | **0.0305 %** | **0.0696 %** | **0.0230 %** |
| FFT, µmol/L | 7.127 | 6.639 | 2.849 |
| FFT, **mol% of ARP** | **0.0356 %** | **0.0332 %** | **0.0142 %** |
| 2-furfural (ARP alone), mol% | 0.0836 % | 0.0697 % | 0.0227 % |
| 2-acetylthiazole (ARP+Cys), mol% | 0.000096 % | 0.00046 % | **0.0229 %** |

**MFT / FFT mass ratio**: pH 6 → **0.857**; pH 7 → **2.096**; pH 8 → **1.616**.
(Both compounds share MW 114.17, so this is also the molar ratio.)
Cross-check against Hofmann & Schieberle 1998 Table 1 (ribose/cys, pH 5.0, 145 °C, SIDA):
MFT/FFT = 19.8/12.1 = **1.636**. **Same order, same sign (MFT ≥ FFT) at pH 7–8** — a real,
independent, cross-lab, cross-method agreement on the branch ratio. **At pH 6 it inverts.**

**★ DIMERISATION SINK FRACTION (molar, derived)** — how much thiol is tied up as its own disulfide:

| pH | disulfide µmol/L | monomer µmol/L | mol dimer / mol monomer | **thiol-equivalents in dimer / free monomer** (×2) |
|---|---:|---:|---:|---:|
| 6 | bis(2-Me-3-furyl)disulf. 0.2638 | MFT 6.105 | 0.0432 | **8.6 %** |
| 7 | 0.4533 | MFT 13.914 | 0.0326 | **6.5 %** |
| 8 | 0.2212 | MFT 4.604 | 0.0480 | **9.6 %** |
| 6 | bis(2-furfuryl)disulf. 0.1561 | FFT 7.127 | 0.0219 | **4.4 %** |
| 7 | 0.0781 | FFT 6.639 | 0.0118 | **2.4 %** |
| 8 | ND | FFT 2.849 | 0 | **0 %** (below LOD) |
| 6 | di-2-furfuryl sulfide 0.0337 | FFT 7.127 | 0.0047 | 0.95 % |
| 7 | 0.0447 | FFT 6.639 | 0.0067 | 1.35 % |
| 8 | ND | FFT 2.849 | 0 | 0 % |

**This is a second, independent measurement of the thiol→disulfide sink that Zhang 2024 measured
at 115 °C in the VB1/Xyl system.** Here it is 6.5–9.6 % of MFT (thiol-equivalents) across
pH 6–8 at 120 °C / 60 min, and it is **near-invariant in pH** (spread 6.5–9.6 %, i.e. ±20 % of
its own mean) while MFT itself swings **3.0×** across the same pH range. That is a strong
structural statement: **the dimerisation branch ratio is roughly pH-independent; the MFT level
is not.** The paper says the same qualitatively: "the variations of their contents with pH were
consistent with their precursors under different pH conditions" (p. 2476).

---

## 3. FIGURE 2 (p. 2476) — FINAL pH. THE BUFFER-FREE WARNING, QUANTIFIED. **[fig, 300 dpi]**

Caption verbatim: "**Final pH value of thermally treated ARP with or without Cys at 120 °C for
60 min under different initial pH values (6, 7, and 8) (the concentrations of ARP and Cys were
both 20 mmol/L).**"

| system | initial pH 6 | initial pH 7 | initial pH 8 |
|---|---:|---:|---:|
| **ARP alone** | **3.99** ᵇ | **4.13** ᵇ | **4.99** ᶜ |
| **ARP + Cys** | **3.22** ᵃ | **3.42** ᵃ | **5.07** ᶜ |

**[fig] read-off, y-axis 0–6 with ticks at 0/2/4/6; my digitisation uncertainty ≈ ±0.06 pH unit.**
Duncan letters a/b/c are printed above the bars and are transcribed exactly.

**★ THIS IS THE MOST IMPORTANT METHODOLOGICAL NUMBER IN THE PAPER.** The "pH 6 / 7 / 8" labels
on Table 1 are **initial** pH of an **unbuffered** system. By the end of the 60-min hold:

- pH 6 → 3.99 (ARP), 3.22 (ARP+Cys): **a drop of 2.0 / 2.8 units.**
- pH 7 → 4.13, 3.42: **a drop of 2.9 / 3.6 units.**
- pH 8 → 4.99, 5.07: **a drop of 3.0 / 2.9 units.**

The pH 6 and pH 7 runs **converge to within 0.2 units of each other** by the end. Whatever
separates the pH 6 and pH 7 columns of Table 1 was decided **early in the ramp**, not at the
endpoint pH. Any model that reads Table 1 as a steady-state pH response will be fitting a
trajectory-integrated quantity to a nominal label.

⚠️ **Consequence for the model.** The Zhou pH series is NOT comparable to Hofmann & Schieberle
1998's pH series, which used **0.5 mol/L phosphate buffer** and therefore held pH. Zhou's is a
**pH-trajectory** experiment; Hofmann's is a **pH-clamped** experiment. Recording them as the
same kind of constraint would be an error. **If the model has no pH-trajectory state, Zhou's
Table 1 must be ingested as an ORDINAL/DIRECTIONAL constraint, not an absolute-pH one.**

The paper's own explanation (p. 2476): "pH decreased after heat treatment, which was mainly
because of the depletion of the ARP and the formation of organic acids from the oxidation of
aldehydes." And the notable exception, verbatim: "when the initial pH was 7, after co-heating of
ARP and Cys, the final pH detected was **significantly lower** compared with that of ARP heated
alone while the content of pyrazines was **higher** with added Cys"; and "When the initial pH was
adjusted to 8... the final pH of ARP and Cys co-heated was **not lower** than that of the ARP
heated alone (Figure 2), which was more favorable for the generation of amino precursors."

---

## 4. ★ FIGURE 3 (p. 2477) — THE Cys+MGO TIME × pH SERIES. **[fig, 300 dpi]**

Caption verbatim: "**Main volatile flavor compounds formed from thermally treated Cys and MGO at
different initial pH values (6, 7, and 8) for 0−90 min (the reaction temperature was 120 °C).**"

All 8 panels are µg/L on the y-axis, reaction time (min) on the x. Series: pH 6 (green circles),
pH 7 (red squares), pH 8 (blue triangles). Error bars = SD, n = 3. **t = 0 is drawn at 0 for
every series in every panel; it is the plotted origin, not an independent measurement.**

Digitisation uncertainty stated per panel as a fraction of the axis full scale; read from
`k3img/z_fig3_full.png` (panels a–d) and `k3img/z_fig3_bot.png` (panels e–h).

### (a) 2,5-dimethylpyrazine — y axis 0–250, ticks every 50. ±3 µg/L
| min | pH 6 | pH 7 | pH 8 |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 |
| 30 | 43 | 59 | 103 |
| 60 | 67 | 94 | 205 |
| 90 | 66 | 106 | **237** |

### (b) 2,6-dimethylpyrazine — y axis 0–30, ticks every 6. ±0.4 µg/L
| min | pH 6 | pH 7 | pH 8 |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 |
| 30 | 5.2 | 7.4 | 10.8 |
| 60 | 7.0 | 10.5 | 24.6 |
| 90 | 7.9 | 12.6 | **25.9** |

⚠️ Paper's own reading of (a) vs (b), verbatim: "Since the aldehyde group of MGO was more
reactive than the ketone group, the content of 2-aminopropanal, which further participated in
the condensation to form 2,6-dimethylpyrazine was less than that of 1-aminopropan-2-one, thereby
the content of 2,5-dimethylpyrazine was higher than that of 2,6-dimethylpyrazine."
**Measured 2,5-DMP / 2,6-DMP ratio at 90 min: pH 6 → 8.4; pH 7 → 8.4; pH 8 → 9.2.**
A clean, pH-insensitive **regiochemical branch ratio of ≈ 8.4–9.2 : 1** for the two dimethyl-
pyrazine isomers out of a *single* α-dicarbonyl. This is directly ingestable as a
**stoichiometric/branching constraint** and needs no absolute-yield calibration to be usable.

### (c) thiophene — y axis 0–20, ticks every 5. ±0.4 µg/L
| min | pH 6 | pH 7 | pH 8 |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 |
| 30 | **0 (ND)** | 0 (ND) | 0 (ND) |
| 60 | 11.4 | 0 (ND) | 0 (ND) |
| 90 | 17.6 | 0 (ND) | 0 (ND) |

**Thiophene is detected ONLY at pH 6, and only from 60 min onward.** A hard on/off in pH — the
sharpest pH switch in the paper. (The pH 8 series is drawn on the zero line and is visually
coincident with pH 7; both read as flat zero.) ⚠️ **There is an induction period**: nothing at
30 min, 11.4 µg/L at 60 min. Zero-order extrapolation through the origin would be wrong here.

### (d) 2-methylthiophene — y axis 0–20, ticks every 5. ±0.4 µg/L
| min | pH 6 | pH 7 | pH 8 |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 |
| 30 | 5.7 | 4.9 | 2.7 |
| 60 | 13.5 | 5.7 | 5.2 |
| 90 | 16.9 | 12.8 | 7.5 |

Monotone **decreasing in pH** at every time point. ⚠️ The pH 7 curve has a visible plateau
between 30 and 60 min then resumes — non-monotone in *shape*, not in level.

### (e) 2-acetylthiazole — y axis 0–800, ticks every 200. ±10 µg/L
| min | pH 6 | pH 7 | pH 8 |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 |
| 30 | 238 | 310 | 425 |
| 60 | 285 | 350 | 665 |
| 90 | 292 | 458 | **755** |

**Cross-system consistency check (mine):** System A (ARP + Cys, pH 8, 120 °C, **60 min**) gave
2-acetylthiazole = **582.34 µg/L** (Table 1). System B (Cys + MGO, pH 8, 120 °C, **60 min**) gives
**665 µg/L [fig]**. **Agreement to 14 %** across two entirely different feedstocks. That is
strong evidence that **MGO is the rate-controlling intermediate for 2-acetylthiazole** and that
the ARP's only role is to supply it — exactly the paper's thesis. **This is a genuine step-level
validation opportunity**: a model that routes ARP → MGO → 2-acetylthiazole must reproduce both
numbers with the same downstream parameters.

### (f) 2-methylthiazole — y axis 0–4, ticks every 1. ±0.06 µg/L
| min | pH 6 | pH 7 | pH 8 |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 |
| 30 | 0.15 | 0.05 | 0.02 |
| 60 | 2.23 | 0.92 | 0.50 |
| 90 | **2.82** | 1.35 | 0.66 |

**Decreasing in pH — the OPPOSITE of 2-acetylthiazole (e), in the same reaction, from the same
two reagents.** ⚠️ And note the same **induction period** as thiophene: almost nothing by 30 min,
then a jump. Paper's mechanism (p. 2478): "nitrogen analogue cysteamine could react with
acetaldehyde to yield 2-methythiazole, and the content of 2-methythiazole at lower pH was higher
than that under alkaline conditions (Figure 3f)."

⚠️ **Cross-system CONFLICT for 2-methylthiazole.** Table 1 (System A, ARP+Cys, 60 min) gives
0.03 → 0.21 → 3.23 µg/L for pH 6 → 7 → 8: **INCREASING in pH.** Figure 3f (System B, Cys+MGO,
60 min) gives 2.23 → 0.92 → 0.50: **DECREASING in pH.** Same compound, same temperature, same
duration, same pH labels, **opposite sign**. The paper does not reconcile this. Do **not** ingest
a 2-methylthiazole pH direction from either system alone.

### (g) pyrazine — y axis 0–50, ticks every 10. ±0.8 µg/L
| min | pH 6 | pH 7 | pH 8 |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 |
| 30 | 37.0 | 20.7 | 12.7 |
| 60 | 46.3 | 19.4 | 14.4 |
| 90 | **47.0** | 20.8 | 15.4 |

### (h) methylpyrazine — y axis 0–200, ticks every 50. ±3 µg/L
| min | pH 6 | pH 7 | pH 8 |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 |
| 30 | 95 | 66 | 65 |
| 60 | 146 | 86 | 106 |
| 90 | **160** | 97 | 117 |

⚠️ **(h) is NON-MONOTONE IN pH**: at 60 and 90 min, pH 8 > pH 7. The pH 6 series is highest
throughout, but the pH 7 and pH 8 series **cross** between 30 and 60 min.

### ★★ THE PARADOX THE PAPER ITSELF FLAGS AND CANNOT EXPLAIN

**In the SAME Cys+MGO reaction: the two DIMETHYL-pyrazines (a, b) INCREASE with pH, while the
unsubstituted pyrazine and methylpyrazine (g, h) DECREASE with pH.** Verbatim, p. 2478:

> "Surprisingly, pyrazine and methylpyrazine were detected in the reaction between Cys and MGO,
> and their formation was favored at lower pH in general (Figure 3g,h). GO was an important
> precursor for the formation of these compounds. It had been reported that GO and MGO could be
> converted into each other in the presence of ammonium sulfide. Therefore, it was reasonable to
> speculate that the reducing activity of Cys probably promoted the conversion of MGO to GO. As
> we know, Cys could be degraded to H₂S, which could be further used to reduce MGO to
> hydroxyacetone. Subsequently, the aldol condensation of MGO to hydroxyacetone led to the
> formation of 2,5-dioxo-3,4-dihydroxyhexane, and this α-dicarbonyl compound might undergo retro
> aldol cleavage at different sites to produce GO (Figure 4). When the pH was low, the Strecker
> degradation became weaker, as well as the reducing activity of H₂S was enhanced, which could
> explain why the content of pyrazine and methylpyrazine was higher at low pH. **Further studies
> are needed to explain this phenomenon.**"

**Structural consequence: "pyrazines" cannot be modelled as one family with one pH term.**
The C0/C1 pyrazines and the C2 pyrazines have **opposite pH signs in the same pot**. Any
family-level pyrazine pH coefficient is falsified by this single figure. This is the pyrazine
analogue of Hofmann & Schieberle 1998's pentose/hexose sign-crossing.

---

## 5. ★ THE H₂S-AVAILABILITY MECHANISM — AS THEY STATE IT, VERBATIM

The brief asked specifically for this. Every sentence in this section is quoted, not paraphrased.

**5.1 The core statement (Abstract, p. 2472):**
> "The formation of thiols and sulfides or 2-acetylthiazole and pyrazines induced by Cys during
> thermal degradation of ARP was pH-dependent. **At low pH levels, the hydrolysis of Cys to
> hydrogen sulfide (H₂S) was promoted**, giving rise to the increase of thiols and sulfides with
> an obvious meaty aroma. However, **alkaline conditions were beneficial for enhancing the
> cyclization or transformation of imine to the enol structure**, which strengthened the formation
> of 2-acetylthiazole and pyrazines with a roasted and nutty aroma. **The imine was derived from
> the nucleophilic addition of Cys and methylglyoxal (MGO) and subsequent decarboxylation.**"

**5.2 The proton mechanism (p. 2478):**
> "H₂S, as an important precursor for the formation of thiophenes and thiols, could be generated
> through the degradation of Cys in the Maillard reaction or through hydrolysis. **Under acidic
> conditions, the thiol group can easily obtain H⁺, and thus the release of H₂S was favored at
> lower pH.** Consequently, the content of thiophene and 2-methylthiophene detected in the system
> between Cys and MGO declined with increasing pH (Figure 3c,d)."

**5.3 The α-dicarbonyl-catalysed hydrolysis route (p. 2478) — this is the mechanistically
specific claim, and it is a *catalysed* route, not spontaneous hydrolysis:**
> "Apart from the Strecker degradation of Cys, **α-dicarbonyl compounds could also become the
> catalysts for the hydrolysis of Cys**. Under acidic conditions, **the protonation of the imine in
> Compound 1 enhanced the reactivity toward H₂O, resulting in the hydrolysis of the Schiff
> structure** to form the corresponding aldehydes and amines (Figure 4)."

**5.4 The competing alkaline branch (p. 2478):**
> "the key regulatory point of 2-acetylthiazole and pyrazine formation was **Compound 1 formed by
> the nucleophilic addition of Cys and MGO followed by decarboxylation** (Figure 4). **The alkaline
> conditions were beneficial for improving the basicity of the imine in Compound 1, thereby
> increasing the polarization degree of π-electrons of the Schiff base. Subsequently, the
> nucleophilic attack of the thiol group on the Schiff base was promoted to form a ring, leading
> to the increase of 2-acetyl-2-thiazoline which was then oxidized to 2-acetylthiazole.**
> Meanwhile, the alkaline environment could also promote the increase of the electron density of
> carbonyl oxygen in the structure of Compound 1, which in turn improved its ability to bind
> protons, as well as promoted the movement of equilibrium toward the enol structure (Figure 4).
> Consequently, the formation of pyrazines was facilitated."

**5.5 The MFT rise-then-fall, explained by H₂S availability (p. 2476) — THE key passage:**
> "Another important thiol produced in the system was 2-methyl-3-furanthiol, which had been
> confirmed to be formed from H₂S and 4-hydroxy-5-methyl-3(2H)-furanone derived from
> 1-deoxypentosone in the Maillard reaction, **and its content increased followed by a decrease
> with increasing pH (Table 1). This could be explained by that the lower acidic conditions
> favored the 1,2-enolization of ARP, resulting in formation of a less amount of
> 1-deoxypentosone, thereby the content of 2-methyl-3-furanthiol was higher at pH 7 in comparison
> to pH 6. When the pH increased to 8, although the formation of 1-deoxypentosone was enhanced,
> the availability of H₂S was inhibited, resulting in a subsequent decline of
> 2-methyl-3-furanthiol.**"

**★ THIS IS A TWO-FACTOR PRODUCT LAW, STATED EXPLICITLY BY THE AUTHORS:**

> **[MFT] ∝ [1-deoxypentosone](pH) × [H₂S available](pH)**
> where **[1-deoxypentosone] increases with pH** (2,3-enolisation route favoured at higher pH)
> and **[H₂S] decreases with pH** (protonation-driven release favoured at low pH).
> The product of a rising and a falling factor gives a **maximum**, and the measured maximum is
> at **pH 7**.

This is *exactly* the kind of structural constraint a mechanistic model should reproduce **for
free**, without a fitted pH term. If the repo's network carries both an enolisation-branch pH
dependence and an H₂S-release pH dependence, **the pH-7 MFT maximum should fall out.** If it
does not, one of the two factors is missing or has the wrong sign. **Treat this as a
first-class structural test, not as three benchmark rows.**

**5.6 The corroborating FFT statement (p. 2476)** — FFT has only ONE pH factor and is therefore
monotone:
> "Owing to the fact that acidic condition was conducive to the production of 2-furfural, the
> synthesis of 2-furfurylthiol decreased with the increase of pH (Table 1)."

Confirmed by the data: 2-furfural (in ARP alone) falls 1606.92 → 1339.37 → 436.63 µg/L across
pH 6→7→8, and FFT falls 813.65 → 757.97 → 325.22. **FFT tracks its furfural precursor
monotonically; MFT does not track its furanone precursor monotonically, because H₂S gates it.**
**The MFT/FFT contrast within a single experiment is the single cleanest piece of evidence in
the corpus that H₂S availability is a distinct, pH-dependent, rate-limiting state variable.**

**5.7 The 2-furfural + H₂S consumption statement (p. 2476)** — furfural is consumed, not just
formed:
> "during the heating of ARP with the addition of Cys, the contents of furans were reduced,
> especially 2-furfural. Although 2-furfural had little contribution to meat flavor, it was a key
> intermediate that could further react with **H₂S provided from Cys** at high temperatures to
> form other flavor compounds such as 2-furfurylthiol."

**Quantified**: adding Cys drops 2-furfural by **62 %** at pH 6 (1606.92 → 608.41), **65 %** at
pH 7 (1339.37 → 473.07), **48 %** at pH 8 (436.63 → 224.98). ⚠️ **Mass-balance check (mine):**
the 2-furfural *lost* at pH 6 is 998.51 µg/L = 10.39 µmol/L; the FFT *gained* is 7.13 µmol/L.
**69 % of the lost furfural is accounted for by FFT alone.** At pH 7: lost 866.30 µg/L =
9.02 µmol/L, FFT gained 6.64 µmol/L → **74 %**. At pH 8: lost 211.65 µg/L = 2.20 µmol/L, FFT
gained 2.85 µmol/L → **129 %, i.e. over-accounted**. The pH 8 point does not close; either
another furfural source opens at pH 8 or the SPME response drifts. **Flag: the furfural→FFT
carbon balance closes at ~70 % at pH 6–7 and fails at pH 8.**

**5.8 The 2,3-butanedione suppression (Table 1) — a Cys competition constraint.** 2,3-Butanedione
in ARP alone rises 8.93 → 33.63 → 126.18 µg/L with pH; adding Cys drives it to **ND, ND, 20.00**.
Complete suppression at pH 6 and 7, **84 % suppression** at pH 8. Cys out-competes whatever forms
2,3-butanedione for the α-dicarbonyl pool. Cross-reference: the authors' own previous paper
(ref 10, Zhou et al. JAFC 2022, 70, 15202−15212, "Competitive formation of 2,3-butanedione and
pyrazines through intervention of added cysteine during thermal processing of alanine-xylose
Amadori compound") is the dedicated study of exactly this. **[RETRIEVE]**

**5.9 The pyrazine competition statement (p. 2476):**
> "the content of pyrazines in the presence of Cys was higher than that formed from the ARP
> heated alone (Table 1), which was due to the fact that **Cys could compete with the regenerated
> alanine to capture α-dicarbonyl compounds** leading to the increase of pyrazines. In addition,
> **Cys could undergo thermal degradation to release ammonia**, which subsequently could react with
> some reactive intermediates such as glycolaldehyde generated from the degradation of ARP to form
> more α-aminoketones, and thus promote the formation of pyrazines."

**5.10 Figure 4 (p. 2478)** is the proposed Cys–MGO pathway scheme (structures only, no numbers).
Nodes named in the text: Cys + MGO → **nucleophilic addition** → **decarboxylation** →
**Compound 1 (the imine / Schiff base)**; then either (acid) **protonation → hydrolysis →
aldehydes + amines + H₂S**, or (base) **thiol attack on Schiff base → ring closure →
2-acetyl-2-thiazoline → [oxidation] → 2-acetylthiazole**, or (base) **enol tautomer → pyrazines**.
Side route: **H₂S reduces MGO → hydroxyacetone; MGO + hydroxyacetone aldol →
2,5-dioxo-3,4-dihydroxyhexane → retro-aldol → GO → pyrazine/methylpyrazine.**
Also named: **mercaptoacetaldehyde** = the Strecker aldehyde of Cys → thiophenes;
**MGO + acetaldehyde aldol → + H₂S → dehydration/reduction → 2-methylthiophene** (ref 33).

---

## 6. ★★ THE CONFLICT WITH HOFMANN & SCHIEBERLE 1998 — DO NOT SUPPRESS IT

| | Hofmann & Schieberle 1998 | Zhou et al. 2023 |
|---|---|---|
| feedstock | free ribose (100 mM) + cysteine (33 mM), **1:3 cys:sugar** | **purified Amadori** (20 mM) + cysteine (20 mM), **1:1** |
| medium | **0.5 mol/L phosphate buffer** (pH held) | **deionized water, unbuffered** (pH falls 2–3.6 units) |
| T / t | **145 °C / 20 min** | **120 °C / 60 min** |
| quantification | **SIDA** (²H-labelled internal standards) | **HS-SPME + external calibration**, 1,2-DCB internal std |
| pH range | 3.0, 5.0, 7.0 | 6, 7, 8 (initial) |
| **MFT vs pH** | **monotone DECREASING with pH**: 55.3 → 19.8 → 2.5 µg/100 mL over pH 3→5→7 | **NON-MONOTONE, MAXIMUM AT pH 7**: 697 → **1589** → 526 µg/L over pH 6→7→8 |
| MFT at the shared pH-7 point | 2.5 µg/100 mL = **25 µg/L** | **1589 µg/L** |

**The overlapping region is pH 6–7 and the two disagree on the SIGN of d[MFT]/dpH.** Hofmann has
MFT falling from pH 5 to pH 7 (19.8 → 2.5, a factor of 7.9 down). Zhou has MFT rising from pH 6
to pH 7 (697 → 1589, a factor of 2.3 up). The absolute levels differ by **64×** at pH 7.

**Candidate reconciliations, none of which the corpus can currently decide between:**
1. **Feedstock state.** Hofmann's system must first *form* the Amadori/deoxyosone from free
   ribose; Zhou's is fed it. If the pH-limiting step in Hofmann is the *condensation* (favoured
   at low pH for the 1,2-enolisation route) and in Zhou it is the *fragmentation*, the two curves
   are measuring different steps and both signs can be right.
2. **pH trajectory.** Zhou's pH-7 run ends at final pH **3.42** — i.e. it spends most of its time
   in *Hofmann's acidic regime*, while Zhou's pH-6 run ends at **3.22**, only 0.2 units lower.
   The 2.3× MFT difference between Zhou's "pH 6" and "pH 7" columns is therefore **not** a 1-unit
   pH difference at all; it is a difference in *how long the system spent above pH 4*.
   **This alone could invert the apparent sign.**
3. **Cys:sugar stoichiometry** (1:3 vs 1:1) changes which reagent limits H₂S.
4. **Temperature** (145 vs 120 °C) changes the relative weight of the two competing branches.
5. **Method** (SIDA vs SPME+external curve) can move an absolute level by 10× but should not
   invert a within-paper sign.

**RECOMMENDATION: ingest Zhou's MFT pH shape as a STRUCTURAL/ORDINAL constraint qualified by
"unbuffered, initial pH, Amadori-fed", and Hofmann's as a separate one qualified by "buffered,
clamped pH, free-sugar-fed". Do NOT merge them into one pH response curve. If a single
model must satisfy both, it needs pH trajectory as a state variable — which is exactly the
kind of thing this pair of papers is diagnostic for.**

---

## 7. ★ SI TABLE S2 — THE ODOR-THRESHOLD LIST (the K2 gap, partially filled)

**"Table S2 The OAVs of the volatile compounds in ARP and ARP + Cys systems"**, SI p. 2.
Column header verbatim: **"Threshold in water (µg/L)"**. "–" = not detected.

⚠️ **PROVENANCE: THE SI GIVES NO CITATION FOR ANY THRESHOLD.** There is no reference column, no
footnote, and the main text cites only "very low odor threshold values (Table S2)" (p. 2476)
without a source. **These are CITED-WITHOUT-SOURCE values.** The basis (**"in water"**) *is*
stated, which is more than most papers give, but the numbers themselves are un-anchored.
**Ingest as thresholds-with-declared-basis-but-unknown-provenance; do NOT treat as measured.**

| compound | **threshold in water (µg/L)** | class |
|---|---:|---|
| furan | 4500 | furan |
| furfural (2-furfural) | 3000 | furan |
| 2-acetylfuran | 10000 | furan |
| pyrazine | 75000 | pyrazine |
| methylpyrazine | 60 | pyrazine |
| 2,5-dimethylpyrazine | 80 | pyrazine |
| 2,6-dimethylpyrazine | 400 | pyrazine |
| thiophene | 84 | thiophene |
| 2-thiophenecarboxaldehyde | 5000 | thiophene |
| thiazole | 38 | thiazole |
| 2,4-dimethylthiazole | 18 | thiazole |
| **2-acetylthiazole** | **10** | thiazole |
| **2-methyl-3-furanthiol (MFT)** | **0.005** | thiol |
| **2-furfurylthiol (FFT)** | **0.006** | thiol |
| **bis(2-methyl-3-furyl) disulfide** | **0.00032** | disulfide |

⚠️ Thresholds are **absent** for: 2,3-butanedione, 2-methylthiophene,
5-methyl-2-thiophenecarboxaldehyde, thieno[3,2-b]thiophene, 2-methylthiazole,
di-2-furfuryl sulfide, bis(2-furfuryl) disulfide. Table S2 simply omits those rows.
**15 of 22 compounds have a threshold; 7 do not.**

**Threshold RATIOS (the quantity the model actually needs, and the one that survives basis
differences):**
- **MFT : FFT = 0.005 : 0.006 = 0.833** — MFT is 1.2× more potent than FFT, in water.
- **bis(2-methyl-3-furyl) disulfide : MFT = 0.00032 : 0.005 = 0.064** — the **disulfide dimer is
  15.6× MORE potent than its own monomer**. Combined with §2.1 (only 6.5–9.6 % of MFT-equivalents
  sit in the dimer), the dimer nonetheless carries **more OAV than the monomer**: see Table S2's
  own OAVs, where at pH 7 MFT = 3.18 × 10⁵ and the disulfide = 3.21 × 10⁵. **The "sink" for MFT
  is not an olfactory sink — it is an olfactory amplifier.** This directly qualifies the
  Zhang 2024 dimerisation-sink finding: mass lost to dimer is *not* aroma lost.
- **MFT : 2-acetylthiazole = 0.005 : 10 = 1/2000.**
- **methylpyrazine : pyrazine = 60 : 75000 = 1/1250** — methyl substitution buys 3 orders of
  magnitude of potency in the pyrazine family.

### 7.1 SI Table S2 — the OAVs as printed (dimensionless; = conc./threshold)

| compound | pH 6 ARP | pH 6 ARP+Cys | pH 7 ARP | pH 7 ARP+Cys | pH 8 ARP | pH 8 ARP+Cys |
|---|---:|---:|---:|---:|---:|---:|
| furan | – | – | – | 6.8 × 10⁻³ | 7.8 × 10⁻³ | 0.014 |
| furfural | 0.54 | 0.20 | 0.45 | 0.16 | 0.15 | 0.07 |
| 2-acetylfuran | 3.23 × 10⁻³ | 6.67 × 10⁻⁴ | 6.22 × 10⁻³ | 6.63 × 10⁻⁴ | 9.34 × 10⁻³ | 9.81 × 10⁻⁴ |
| pyrazine | – | – | – | 8.08 × 10⁻⁵ | 5.32 × 10⁻⁵ | 0.0024 |
| methylpyrazine | – | – | – | 0.017 | 0.066 | 1.65 |
| 2,5-dimethylpyrazine | – | – | – | – | 0.0034 | 0.268 |
| 2,6-dimethylpyrazine | – | – | – | – | – | 0.002 |
| thiophene | – | 0.16 | – | 0.26 | – | 0.20 |
| 2-thiophenecarboxaldehyde | – | 4.33 × 10⁻³ | – | 9.12 × 10⁻³ | – | 8.73 × 10⁻³ |
| thiazole | – | 0.14 | – | 0.29 | – | 2.83 |
| 2,4-dimethylthiazole | – | – | – | – | – | 0.043 |
| **2-acetylthiazole** | – | 0.243 | – | 1.17 | – | **58.234** |
| **2-methyl-3-furanthiol** | – | **1.39 × 10⁵** | – | **3.18 × 10⁵** | – | **1.05 × 10⁵** |
| **2-furfurylthiol** | – | **1.36 × 10⁵** | – | **1.26 × 10⁵** | – | **5.42 × 10⁴** |
| **bis(2-methyl-3-furyl) disulfide** | – | **1.87 × 10⁵** | – | **3.21 × 10⁵** | – | **1.56 × 10⁵** |

⚠️ **The OAVs are internally consistent with Table 1 and the thresholds** — I spot-checked four:
MFT pH 7: 1588.57 / 0.005 = 3.177 × 10⁵ ✓ (printed 3.18 × 10⁵).
2-acetylthiazole pH 8: 582.34 / 10 = 58.234 ✓ (printed exactly).
methylpyrazine pH 8 ARP+Cys: 99.00 / 60 = 1.65 ✓.
bis(2-methyl-3-furyl) disulfide pH 6: 59.70 / 0.00032 = 1.866 × 10⁵ ✓ (printed 1.87 × 10⁵).
**No transcription errors detected in the SI.**

**★ The whole aroma is three compounds.** At every pH, MFT + FFT + the MFT disulfide carry
OAVs of 10⁴–10⁵ while everything else is ≤ 58. **The thiol/disulfide triad is 3–7 orders of
magnitude above the rest of the volatile profile.** Any aroma-weighted objective built on this
system is, numerically, a three-compound objective.

## 8. SI TABLE S1 — CALIBRATION EQUATIONS (all 22 compounds)

Form: **y = m·x + c**, where **x = concentration of the flavour compound** and **y = peak-area
ratio (compound / 1,2-dichlorobenzene internal standard)**. The units of x are not restated in
the SI; from Table 1 they must be **µg/L**, so m has units **(L/µg)**.

| compound | calibration equation | R² |
|---|---|---:|
| 2,3-butanedione | y = 0.0109x + 0.0331 | 0.9972 |
| furan | y = 0.0015x + 0.0067 | 0.9951 |
| furfural | y = 0.0191x − 0.3498 | 0.9900 |
| 2-acetylfuran | y = 0.0107x + 0.0903 | 0.9924 |
| pyrazine | y = 0.0012x − 0.0017 | 0.9988 |
| methylpyrazine | y = 0.0035x + 0.0031 | 0.9976 |
| 2,5-dimethylpyrazine | y = 0.0031x + 0.0043 | 0.9981 |
| 2,6-dimethylpyrazine | y = 0.0063x + 0.003 | 0.9963 |
| thiophene | y = 0.0043x + 0.0025 | 0.9996 |
| 2-methylthiophene | y = 0.003x + 0.0038 | 0.9928 |
| 2-thiophenecarboxaldehyde | y = 0.0113x + 0.02 | 0.9959 |
| 5-methyl-2-thiophenecarboxaldehyde | y = 0.0145x + 0.0062 | 0.9982 |
| thieno[3,2-b]thiophene | y = 0.0277x − 0.0029 | 0.9972 |
| 2-methylthiazole | y = 0.0252x + 0.0058 | 0.9942 |
| thiazole | y = 0.0043x + 0.0061 | 0.9933 |
| 2,4-dimethylthiazole | y = 0.0493x + 0.0131 | 0.9967 |
| 2-acetylthiazole | y = 0.0132x + 0.0041 | 0.9975 |
| **2-methyl-3-furanthiol** | **y = 0.001x + 0.0029** | 0.9966 |
| **2-furfurylthiol** | **y = 0.0018x + 0.0037** | 0.9951 |
| difurfuryl sulfide | y = 0.0181x + 0.002 | 0.9947 |
| **bis(2-methyl-3-furyl) disulfide** | **y = 0.0208x − 0.7126** | 0.9984 |
| bis(2-furfuryl) disulfide | y = 0.0051x − 0.021 | 0.9966 |

⚠️ **Three curves have NEGATIVE intercepts**: furfural (−0.3498), bis(2-methyl-3-furyl)
disulfide (−0.7126), bis(2-furfuryl) disulfide (−0.021), pyrazine (−0.0017),
thieno[3,2-b]thiophene (−0.0029). For the disulfide the intercept is large relative to the
slope: y = 0 at x = 34.3 µg/L, i.e. **the calibration curve predicts zero signal below ~34 µg/L
of bis(2-methyl-3-furyl) disulfide**. The pH 8 value in Table 1 is **50.07 µg/L** — only 1.5×
above that pseudo-LOD. ⚠️ **The disulfide numbers, and hence my derived dimerisation-sink
fractions in §2.1, carry a systematic uncertainty that the stated ±SD does not capture.**
Similarly bis(2-furfuryl) disulfide zeroes at x = 4.1 µg/L, and "ND" at pH 8 may mean
"below 4 µg/L", not "absent".

⚠️ **No concentration RANGE, no number of calibration levels, and no LOD/LOQ is given for any
curve.** Extrapolation beyond the (unstated) calibrated range cannot be ruled out for the
largest values — e.g. MFT at 1588.57 µg/L.

---

## 9. THE ELECTRONIC-NOSE / PCA RESULT (Figure 1, p. 2474) — qualitative, one number

> "PC1 explained **99.86 %** of the total variation, and **0.10 %** of the total variation was
> contributed by PC2. The sum of the two variance contribution rates accounted for **99.96 %**."

> "the samples were divided into **6 clusters**... The non-Cys-added group and the Cys-added group
> both showed **negative progress along the abscissa axis (PC1) as the pH rose**. Furthermore, at
> **pH 6 and 7**, the groups of ARP heated alone and the groups of heated ARP with Cys involvement
> were **clearly separated**... However, while the pH was adjusted to **8**, no matter Cys addition
> or not, the groups of thermal processed ARP were **both in the third quadrant at positions
> closer to each other, indicating that there was no significant distinction** between the flavor
> profiles of the two."

Heracles II fast GC e-nose, two columns (MXT-5 and MXT-1701, both 10 m × 0.18 mm), two FIDs;
headspace 60 °C / 30 min; 5000 µL injected; ramp 50 °C (2 s) → 80 °C at 1 °C/s → 250 °C at
2 °C/s, hold 60 s. n = 4.

**Directional claim only: "Cys addition changes the flavour profile at pH 6 and 7, but NOT at
pH 8."** ⚠️ A 99.86 % PC1 on FID peak areas is dominated by total signal magnitude; treat the
clustering as weak evidence.

---

## 10. WHAT THIS PAPER UNLOCKS

1. **A fed-Amadori initial condition.** The model can now be initialised at a *named, purified,
   structurally characterised* Amadori compound at a known 20 mM, with the sugar removed by ion
   exchange. This is the first such initial condition in the corpus for the pentose/alanine lane.
   Every prior "Amadori" constraint has been inferred from a sugar + amino acid mixture.
2. **A step-level MGO node validation.** The 2-acetylthiazole agreement between System A
   (582 µg/L, ARP-fed) and System B (665 µg/L, MGO-fed) at the same pH/T/t constrains the
   ARP → MGO flux to within ~15 % **if** the downstream MGO → 2-acetylthiazole step is shared.
   That is a real, checkable, two-system constraint on one node.
3. **A mechanistic, author-stated product law for MFT vs pH** (§5.5) that a correct network
   should reproduce without a fitted pH term — and a pH-7 maximum to test it against.
4. **A pH-invariant thiol dimerisation branch ratio** (6.5–9.6 % of MFT-equivalents; §2.1),
   independently corroborating Zhang 2024's MFT→MFT-MFT sink at a different temperature,
   different feedstock, different lab.
5. **A falsifier for family-level pyrazine pH terms** (§4): C0/C1 pyrazines and C2 pyrazines have
   opposite pH signs in the same pot.
6. **A quantified pH-drift warning** (§3) that retro-actively qualifies *every* unbuffered
   Maillard pH series in the corpus.
7. **A 15-compound threshold list with a declared basis** (§7), and the ratio
   bis(2-methyl-3-furyl)disulfide : MFT = 0.064 that turns the "MFT dimerisation sink" from an
   aroma loss into an aroma gain.
8. **Two [RETRIEVE] pointers**, both by the same group, both directly on-topic:
   - **ref 10** — Zhou, Xia, Cui, Hayat, Zhang, Ho, *JAFC* **2022**, 70, 15202−15212,
     "Competitive formation of 2,3-butanedione and pyrazines through intervention of added
     cysteine during thermal processing of alanine-xylose Amadori compound". **The dedicated
     study of the Cys/α-dicarbonyl competition that this paper only measures in passing.**
   - **ref 14** — Liu, Yu, Zhou, Xu, Hayat, Zhang, Ho, *JAFC* **2022**, 70, 11643−11651,
     "Formation priority of pyrazines and 2-acetylthiazole dependent on the added cysteine and
     fragments of deoxyosones during the thermal process of the glycine-ribose Amadori compound".
     **The MGO-vs-GO selectivity study — the missing half of §4's paradox.**

---

## 11. INGESTION RECOMMENDATION

| item | class | note |
|---|---|---|
| **Table 1, all 22 compounds × 6 conditions, µg/L** | **benchmark set (absolute, qualified)** | HS-SPME + external curve, **not** SIDA; **unbuffered, initial pH**; do not mix with SIDA rows in the same residual pool |
| **MFT pH-7 maximum (697 / 1589 / 526)** | **★ STRUCTURAL constraint, ordinal** | must emerge from the two-factor law in §5.5; **do not fit a pH term to reproduce it** |
| **FFT monotone decline with pH (814 / 758 / 325)** | **directional, high confidence** | tracks its 2-furfural precursor; the control case against MFT |
| **MFT/FFT ratio 0.86 / 2.10 / 1.62** | **benchmark (branch ratio)** | ratio within one run; method-independent; cross-validates against Hofmann 1998's 1.64 at pH 5 |
| **Dimer/monomer 6.5–9.6 mol% (MFT), 0–4.4 % (FFT)** | **benchmark (sink fraction)** | ⚠️ qualified by the disulfide calibration intercept (§8): pseudo-LOD ≈ 34 µg/L |
| **2,5-DMP / 2,6-DMP = 8.4–9.2** | **benchmark (regiochemical branch ratio)** | pH-insensitive; from a single α-dicarbonyl; very clean |
| **Figure 3, 8 compounds × 3 pH × 3 times, µg/L** | **benchmark set (time course, fed-MGO)** | **[fig]** digitised, uncertainties stated per panel; t = 0 is a plotted origin not a datum |
| **2-acetylthiazole 582 (System A) vs 665 (System B), pH 8, 60 min** | **★ step-level cross-system constraint** | the single best node test in the paper |
| **Figure 2 final pH (3.99/4.13/4.99; 3.22/3.42/5.07)** | **benchmark (pH trajectory endpoint)** | ⚠️ **must be ingested alongside Table 1 or Table 1 is misleading** |
| **Table S2 thresholds, 15 compounds, µg/L in water** | **threshold list, basis declared, PROVENANCE UNKNOWN** | no citation anywhere in paper or SI; ingest as cited-without-source |
| **Threshold ratios (MFT:FFT 0.83; disulfide:MFT 0.064)** | **usable, ratio form** | ratios are robust to the unknown provenance of the absolute values |
| **Table S1 calibration equations** | **method metadata** | ingest for the LOD/intercept warnings, not as chemistry |
| **Table S2 OAVs** | **derived, verified** | = Table 1 / Table S2 threshold; I re-checked 4 of 15, all exact |
| **2-methylthiazole pH direction** | **DO NOT INGEST — systems disagree in sign** | Table 1 increasing, Fig. 3f decreasing, same T/t/pH labels |
| Absolute µg/L compared against Hofmann 1998 SIDA values | **DO NOT** | 64× apart at pH 7; different method, feedstock, buffer, T |
| Any rate constant, Ea, or Arrhenius parameter | **NOT PRESENT — single temperature (120 °C)** | Fig. 3 is a time series but the paper fits nothing; any k is *your* fit |
| PCA / e-nose clusters | **weak directional only** | PC1 = 99.86 % is a magnitude axis |
| **ref 10** (JAFC 2022, 70, 15202) and **ref 14** (JAFC 2022, 70, 11643) | **[RETRIEVE]** | same group; the dedicated competition and selectivity studies |
