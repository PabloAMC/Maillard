# Meynier et al. 2002/2003 — COMPLETE TRANSCRIPTION
### Air/matrix partition of hexanal, t-2-hexenal and three esters over water, skim milk, anhydrous milk fat and full-fat cream, 30–80 °C.

**Full extraction of every number in `data/articles/Meynier2002.pdf`.**
Wave K4b, 2026-08-28. Read-only extraction: no repo file outside
`data/lit/extraction_dossiers/` was written or modified.

Table III and Table IV re-read from a **400 dpi render**; the text layer is clean but the
`×` and superscript characters are mangled by the 2003 Distiller (`4.4 ´ 10–4` for
`4.4 × 10⁻⁴`), so every exponent below was confirmed against the raster.

---

## 0. PAPER IDENTITY — ⚠️ FILENAME YEAR IS WRONG BY ONE. Content is the expected paper.

| field | value as printed |
|---|---|
| Authors | **Anne Meynier, Aurélie Garillon, Laurent Lethuaut, Claude Genot** |
| Title | **"Partition of five aroma compounds between air and skim milk, anhydrous milk fat or full-fat cream"** |
| Venue | ***Lait* 83 (2003) 223–235** — © INRA, EDP Sciences |
| DOI | **10.1051/lait:2003012** |
| Dates | **Received 17 September 2002; accepted 11 December 2002**; published in the 2003 volume |
| Affiliation | **INRA, LEIMA, BP 71627, 44316 Nantes Cedex 3, France** |
| Funding | Région des Pays de la Loire, VANAM programme *"Interfaces, Emulsions et Mousses Liquides"* |
| **⚠️ file naming** | The repo file is `Meynier2002.pdf`. **The paper's publication year is 2003** (*Lait* 83, 2003). 2002 is the submission/acceptance year. **Cite as Meynier et al. 2003 or the DOI will not match a `Meynier2002` bibkey.** This is a bibkey hazard, not a wrong-file error — the content is exactly the expected air/matrix partition study. |
| PDF character | 13 pages, clean text layer; **all four data tables fully legible; no cell unreadable.** |

**Provenance codes:** **[M]** measured and printed · **[C]** cited by the authors ·
**[Z]** derived by this wave, never printed · **[F]** fitted by the authors and printed.

---

## 1. THE ONE-PARAGRAPH ANSWER

**Table III (20 partition coefficients at 30 °C, each with an SD) and Table IV (20
enthalpies of affinity, each with an SD) are complete and legible.** But **the paper prints
no partition coefficient at any temperature other than 30 °C** — the 30–80 °C behaviour
lives entirely in Figures 1–4, which are line plots on linear axes with no tabulated values.
**The temperature dependence is nonetheless fully recoverable**, because Table IV's
"enthalpy of affinity" *is* the van 't Hoff slope of exactly those curves, and reconstructing
`K(T) = K(303.15)·exp(−ΔH_aff/R·(1/T − 1/303.15))` **reproduces four of the paper's own
printed prose checkpoints to within 2–8 %** (§5.1). That reconstruction is delivered in §5
as a 5 compound × 4 media × 6 temperature grid. Two findings the paper does not draw:
**(i) its own measured air–water partition coefficients are a geometric-mean 6.2× BELOW the
Henry's-law constants it prints in its own Table I** (§6.1), across all five compounds and
both chemical classes — a systematic calibration offset that makes the *ratios* usable and
the *absolutes* not; and **(ii) t-2-hexenal, the α,β-unsaturated aldehyde, is retained
85 % by skim milk at every temperature while hexanal is retained 28 → 44 %** — a
temperature-independent 2-alkenal penalty that is the partition-side counterpart of the
≈2–3× α,β-unsaturation threshold penalty in `k2_matrix_and_thresholds.md` §D.3.

---

## 2. SYSTEM COMPOSITION — applies to every number in this dossier

Stated on pp. 224–226 (§§2.1–2.4).

| variable | value as printed |
|---|---|
| Compounds (5) | **amyl acetate, isoamyl acetate, ethyl pentanoate, hexanal, t-2-hexenal** |
| Purity | **> 98 %** (all four from Aldrich, Saint-Quentin-Falavier); **isoamyl acetate a gift from IFF, Longvic** |
| **Water** | ultrapure, **Millipore system** |
| **Skim milk (sm)** | **"Lactel" skim milk UHT**, purchased from a local supermarket. **Lipid content stated only as "less than 0.3 %".** Protein content **not measured** — the discussion cites literature values: **caseins 29.5 g/L, β-lactoglobulin 3.2 g/L, α-lactalbumin 1.2 g/L** [C] Swaisgood 1995. pH **not stated**. |
| **Full-fat cream (ffc)** | **"Elle & Vire" full-fat cream, 30 % milk fat**, local supermarket. **Dispersed-phase mass fraction = 0.3.** Droplet size **not measured.** |
| **Anhydrous milk fat (amf)** | **France Beurre (Quimper)**, **melting point 32–35 °C**. Pretreatment: **melted at 60 °C for 1 h to erase thermal history, then held at 40 °C** before flavouring and analysis. |
| Dosing | stock solutions weighed accurately, added **at room temperature**; **five concentrations per aroma per medium**; **all below the maximum solubility in that medium** |
| Equilibration of the dosed liquid | **1 h at room temperature in closed flasks under mild stirring** |
| Headspace vial | **3.000 g of sample in a 22.4 mL headspace vial**, sealed with **PTFE–silicon caps**; **duplicates per solution** |
| **Headspace equilibration** | **30 min at the measurement temperature**, ⚠️ **except skim milk, which required 60 min** (determined experimentally beforehand) |
| Temperatures | **30 to 80 °C** (the figures show **30, 40, 50, 60, 70, 80 °C**) |
| Sampler | **Perkin-Elmer HS 40XL** automatic headspace sampler; **pressurisation 2.0 min**; **injected volume 75 µL to 400 µL depending on temperature and detector sensitivity/linearity** |
| GC | **HP 5890**, splitless, **inlet 250 °C**, purge off for the first **1.5 min**; oven **40 °C for 1 min → 7 °C/min → 140 °C**; column **DB 225, 30 m × 0.32 mm i.d., 0.25 µm film** (J&W); **helium 1.9 mL/min**; **FID at 250 °C** |
| ⚠️ **Calibration** | **"The flame ionisation detector (FID) was calibrated by manual injections of aroma compounds in cyclohexane."** Four concentrations per aroma in duplicate, by three operators; **response-factor CV < 5 % between operators.** **This is a LIQUID-phase calibration used to quantify a GAS-phase sample.** See §6.1 — it is the most likely origin of the systematic 6.2× offset. |
| Replication of each K | **mean of 10 measurements = 5 concentrations × duplicate** |
| Linearity | **"Whatever the temperature, the media, and the aroma compounds, headspace concentration was a linear function of liquid concentration over the studied range."** (p. 226) |
| Statistics | **ANOVA + LSD, P < 0.05 throughout**, Statgraphics Plus 3.0 |
| **Definition of K** | *"the ratio of the concentration of each aroma in the gaseous phase (ng·mL⁻¹) to its concentration in the liquid phase (ng·mL⁻¹)"* — i.e. **dimensionless air/matrix, gas over liquid**, so **larger K = more volatile.** ⚠️ Note the liquid concentration is **"calculated after an accurate weighing"**, i.e. the **nominal** liquid concentration, **not corrected for the amount that left for the headspace.** With 3 g in a 22.4 mL vial (V_HS/V_liq ≈ 6.5) and K ≈ 5×10⁻³, the depletion is ≈3 % — negligible; noted for completeness. |

---

## 3. TABLE I — physicochemical constants. **Anchor: Table I, p. 225.** ALL [C], NONE measured.

Title as printed: *"Physicochemical and thermodynamic characteristics of the five aroma
compounds."* Source footnote as printed: ***`http://esc-plaza.syrres.com/interkow/physdemo.htm`***
(the SRC/Syracuse Research **PHYSPROP / EPI Suite** demo server). Footnote (1): *"calculated
according UNIFAC."* Footnote: *"All data were collected at 25 °C (298 K)."*

| compound | **MW** | **BP °C** | **P_sat mm Hg** | **S g·L⁻¹** | **Log P** | **Henry's constant Atm·m³·mol⁻¹** | **γ (UNIFAC)** |
|---|---:|---:|---:|---:|---:|---:|---:|
| Isoamyl acetate | 130 | 142 | 5.6 | 2.0 | **2.26** | **5.87 × 10⁻⁴** | 3673 |
| Amyl acetate | 130 | 149 | 3.5 | 1.7 | **2.3** | **3.88 × 10⁻⁴** | 3664 |
| Ethyl pentanoate | 130 | 145 | 4.8 | 2.2 | **2.34** | **3.72 × 10⁻⁴** | 5345 |
| **Hexanal** | **100** | **131** | **11.3** | **5.6** | **1.78** | **2.13 × 10⁻⁴** | **82** |
| **t-2-hexenal** | **98** | **146** | **6.6** | **5.3** | **1.58** | **4.89 × 10⁻⁵** | **173** |

⚠️ **Every value in Table I is cited from a public estimation server, none is measured**,
and the γ column is a **UNIFAC calculation**, not an experiment. **Ingest only with
`provenance: cited_from_EPI_Suite_estimator`.** The activity coefficients are the tell:
γ = 3664–5345 for the three esters versus **82 for hexanal** — a 45× gap that is a property
of the UNIFAC group contributions, not of any measurement in this paper.

---

## 4. TABLE II — dosing ranges. **Anchor: Table II, p. 225.** All legible. [M]

Title as printed: *"Ranges of aroma concentrations (mg·L⁻¹) in the different media for
partition coefficient determination."* **Units as printed: mg·L⁻¹.**

| compound | **Water** | **Skim milk** | **Full-fat cream** | **Anhydrous milk fat** |
|---|---|---|---|---|
| Isoamyl acetate | 2–200 | 2–200 | 65–6500 | 200–20000 |
| Amyl acetate | 2–200 | 2–200 | 120–1200 | 100–8000 |
| Ethyl pentanoate | 6–600 | 5–500 | 50–5000 | 150–15000 |
| **Hexanal** | **1–100** | **2–200** | **20–2000** | **40–4000** |
| **t-2-hexenal** | **3–300** | ⚠️ **100–5800** | **20–2000** | **200–20000** |

**⚠️ The t-2-hexenal / skim-milk row is 33–1 933× above the corresponding water row**, and
the paper says why (p. 228): *"the accurate quantification of this compound in the gaseous
phase over skim milk needed much higher concentration in the liquid phase: 5.8 to
0.2 g·L⁻¹ in milk compared with 0.3 to 3 × 10⁻³ g·L⁻¹ in water."* And — critically —
***"In presence of such amounts of t-2-hexenal, the skim milk turned in brownish in colour,
even when the solution was kept for a few hours at 4 °C."***
**⇒ The t-2-hexenal / skim-milk partition coefficient was measured on a browning, visibly
reacting system at up to 5.8 g/L, a concentration 1 933× above the water arm.** The paper
itself concedes the constant is only "apparent". **This is the single largest caveat in the
paper and it must ride on every t-2-hexenal/skim-milk record.**

---

## 5. TABLE III — THE 20 MEASURED PARTITION COEFFICIENTS AT 30 °C. **Anchor: Table III, p. 226.** [M]

Title as printed: *"Measured air-medium partition coefficients (× 10³) of aroma compounds at
30 °C."* Note as printed: *"Values are mean of 10 measurements (5 concentrations in
duplicate); (bold: mean, plain: standard deviation)."*
**All values below are the printed table entries × 10⁻³, i.e. dimensionless air/matrix.**

| medium | **Isoamyl acetate** | **Amyl acetate** | **Ethyl pentanoate** | **Hexanal** | **t-2-hexenal** |
|---|---:|---:|---:|---:|---:|
| **Water** | **4.5 ± 0.6** ×10⁻³ | **3.0 ± 0.4** ×10⁻³ | **2.8 ± 0.1** ×10⁻³ | **2.5 ± 0.1** ×10⁻³ | **0.44 ± 0.07** ×10⁻³ |
| **Skim milk** | **4.2 ± 0.5** ×10⁻³ | **2.5 ± 0.2** ×10⁻³ | **2.1 ± 0.3** ×10⁻³ | **1.8 ± 0.2** ×10⁻³ | **0.064 ± 0.027** ×10⁻³ |
| **Anhydrous milk fat** | **0.057 ± 0.003** ×10⁻³ | **0.014 ± 0.002** ×10⁻³ | **0.013 ± 0.002** ×10⁻³ | **0.082 ± 0.006** ×10⁻³ | **0.010 ± 0.001** ×10⁻³ |
| **Full-fat cream** | **0.14 ± 0.02** ×10⁻³ | **0.060 ± 0.003** ×10⁻³ | **0.037 ± 0.004** ×10⁻³ | **0.16 ± 0.03** ×10⁻³ | **0.018 ± 0.005** ×10⁻³ |

### 5.1 Matrix/water suppression factors at 30 °C, computed by this wave [Z]

| compound | **water/skim milk** | **water/amf** | **water/ffc** | **amf/ffc** |
|---|---:|---:|---:|---:|
| isoamyl acetate | **1.07×** | **78.9×** | 32.1× | 0.41× |
| amyl acetate | **1.20×** | **214×** | 50.0× | 0.23× |
| ethyl pentanoate | **1.33×** | **215×** | 75.7× | 0.35× |
| **hexanal** | **1.39×** | **30.5×** | **15.6×** | 0.51× |
| **t-2-hexenal** | **6.88×** | **44.0×** | **24.4×** | 0.56× |

**Two structural facts fall straight out.**
1. **Skim milk (≈34 g/L total protein, <3 g/L fat) suppresses the headspace by only
   1.07–1.39× for the three esters and hexanal, and by 6.88× for t-2-hexenal.**
   That is a **direct, independent confirmation of `k2_matrix_and_thresholds.md` §B.4's
   ceiling on reversible binding**: at a protein loading (34 g/L) close to Andriot's 4 %
   β-lactoglobulin (40 g/L), the measured headspace suppression is **1.1–1.4×** for
   non-reactive ligands, against Andriot's 1.25–3.7×. Same order, independent lab, real food
   matrix, four compounds. **Reversible protein binding remains a single-digit factor.**
2. **The one compound that breaks the pattern is the α,β-unsaturated aldehyde, at 6.88×,
   and the paper attributes it to COVALENT chemistry, not binding** (p. 228, and again
   pp. 231–232): *"covalent binding with histidyl and lysyl residues of proteins are probably
   involved. This irreversible interaction should not be disrupted during consumption."*
   ⚠️ **The anhydrous milk fat column shows the opposite ranking (t-2-hexenal 44× vs hexanal
   30.5×, a mere 1.4× gap) — lipids do not discriminate the two aldehydes, proteins do by
   4.9×.** That contrast is the cleanest evidence in the paper that the skim-milk t-2-hexenal
   effect is chemistry, not hydrophobicity.

⚠️ **The paper notices, and cannot explain, that the emulsion is a worse solvent than the
pure fat.** For all five compounds `K_ffc > K_amf` (amf/ffc = 0.23–0.56), i.e. **the 30 %-fat
cream holds the aroma LESS tightly than 100 % fat** — as expected — but see §7.3: the
measured cream values also diverge from what the two-phase model predicts, in **both
directions depending on the compound**.

---

## 6. TABLE IV — THE 20 ENTHALPIES OF AFFINITY. **Anchor: Table IV, p. 230.** [F]

Title as printed: *"Enthalpies of affinity of the five aroma compounds in water and various
dairy matrices."* Units as printed: **ΔH (kJ·mol⁻¹)**. Note: *"Bold: mean, plain: standard
deviation."* Definition, eq. (3), p. 230: **`ln K_a,medium = −ΔH_aff/(RT) + C′`**, with
**`ΔH_aff = −(slope × R)`, R = 8.314 J·mol⁻¹.**

| ΔH (kJ·mol⁻¹) | **Isoamyl acetate** | **Amyl acetate** | **Ethyl pentanoate** | **Hexanal** | **t-2-hexenal** |
|---|---:|---:|---:|---:|---:|
| **Pure compound** (ΔH of vaporisation) | **44.0** ⁽ᵃ⁾ | **41.9** ⁽ᵃ⁾ | **46.0** ⁽ᵇ⁾ | **45.3** ⁽ᵃ⁾ | **44.8** ⁽ᶜ⁾ |
| **Water** | **35.7 ± 1.6** | **41.2 ± 1.7** | **40.3 ± 1.3** ⁽¹⁾ | **40.4 ± 1.7** | **47.2 ± 1.8** |
| **Skim milk** | **32.0 ± 0.7** | **36.5 ± 1.7** | **39.7 ± 0.4** ⁽²⁾ | **36.0 ± 1.0** | **49.2 ± 4.5** ⁽²⁾ |
| **Anhydrous milk fat** | **47.2 ± 3.3** | **44.5 ± 1.1** | **55.6 ± 3.3** | **37.9 ± 1.3** | **45.5 ± 3.9** |
| **Full-fat cream** | **41.3 ± 0.9** | **47.7 ± 1.0** | **56.5 ± 1.1** | **36.8 ± 1.0** | **55.3 ± 3.5** |

Footnotes as printed:
⁽ᵃ⁾ *Calculated from data in the Handbook of Chemistry and Physics from the variation of
vapour pressure as a function of temperature* [C] Weast & Astle 1980, 61st edn.
⁽ᵇ⁾ *from Ducros et al. (1980)* [C]. ⁽ᶜ⁾ *from Philippe et al. (2002)* [C].
**⁽¹⁾ calculated from 30 to 50 °C only. ⁽²⁾ calculated from 30 to 60 °C only.**

**⚠️ THE TWO SUPERSCRIPT FOOTNOTES ARE LOAD-BEARING AND EASY TO MISS.** Three of the twenty
ΔH values are **not** fitted over the full 30–80 °C range:
- **ethyl pentanoate / water — 30 to 50 °C only**
- **ethyl pentanoate / skim milk — 30 to 60 °C only**
- **t-2-hexenal / skim milk — 30 to 60 °C only**

The reason is printed on p. 228: *"The partition coefficient over skim milk tended to level
off between 70 °C and 80 °C whereas it increased exponentially for other esters."*
**⇒ For those three series the van 't Hoff form is invalid above the stated ceiling and any
extrapolation to 70–80 °C is wrong.** §7 marks the affected cells.

### 6.1 What the ΔH table says, as the authors read it (pp. 230–231, all their words)

- *"in contrast to the partition coefficients, all data are in the same order of magnitude,
  varying from 32 to 57 kJ·mol⁻¹."*
- *"With the exception of t-2-hexenal and, to a lesser extent, amyl acetate, the so-called
  enthalpies of affinity in water were lower than the enthalpies of vaporisation of pure
  compounds. These findings suggest that aroma-aroma interactions were greater than
  aroma-water interactions."*
- *"Enthalpies of affinity tended to be lower in skim milk than in water, with the exception
  of t-2-hexenal and, to a lesser extent, of ethyl pentanoate. The increase of enthalpy of
  affinity of t-2-hexenal was logical considering the strong interaction of this aroma with
  dairy proteins."*
- *"the differences between the enthalpies of affinity of the esters were by far more
  significant than the differences in log P values."*
- *"The enthalpies of affinity of aldehydes were similar in milk fat and in lipid-free
  media."* (hexanal 37.9 amf vs 40.4 water vs 36.0 sm; t-2-hexenal 45.5 vs 47.2 vs 49.2)
- *"With the exception of isoamyl acetate and hexanal, the enthalpies of affinity in full-fat
  cream gave unexpected results ... Values measured in full-fat cream were greater than or
  similar to those measured over milk, while we expected values intermediate between those in
  anhydrous milk fat and skim milk ... We currently have no clear explanation for this
  phenomenon."*

---

## 7. THE K(T) GRID — RECONSTRUCTED BY THIS WAVE [Z], AND VALIDATED AGAINST THE PAPER'S OWN PROSE

**The paper prints no K above 30 °C.** Combining Table III (the 30 °C anchor) with Table IV
(the van 't Hoff slope) recovers the whole of Figures 1–4:

```
K(T) = K(303.15 K) · exp( −ΔH_aff/R · (1/T − 1/303.15) )
```

### 7.1 The validation — four independent checkpoints the paper prints in prose

| the paper's printed statement | page | **reconstruction [Z]** | agreement |
|---|:--:|---:|---|
| *"Measured partition coefficients over water ranged from 4.4 × 10⁻⁴ for t-2-hexenal at 30 °C to **3.4 × 10⁻² for isoamyl acetate at 80 °C**"* | 227 | **3.34 × 10⁻²** | **✅ 2 %** |
| *"For amyl acetate, retention became significant **above 70 °C (~30 %)**"* | 228 | **29.5 %** at 70 °C | **✅ 2 %** |
| *"the retention [of hexanal] reached **42 % at 80 °C**"* | 228 | **43.8 %** | **✅ 4 %** |
| *"The partition coefficient [of t-2-hexenal] over skim milk was approximately **85 % lower** than over water **whatever the temperature**"* | 228 | **85.5 → 83.7 %**, flat | **✅ and flat** |
| *"retention of isoamyl acetate in skim milk varied from **20 % at 40 °C** to **26 % at 80 °C**"* | 227–228 | **10.9 % → 24.2 %** | **⚠️ 80 °C ✅ (7 %); 40 °C ✗ off by 1.8×** |
| *"[ethyl pentanoate] Retention ranged from **25 % at 30 °C to 48 % at 80 °C**"* | 228 | **25.0 % → 27.5 %** | **⚠️ 30 °C ✅ exact; 80 °C ✗** — but this is the series whose ΔH footnotes ⁽¹⁾⁽²⁾ say the fit stops at 50/60 °C, and whose curve the paper says "levels off". **The reconstruction is invalid there by the paper's own footnote, not by error.** |

**⇒ The reconstruction is validated for 4 of the 5 water series and 3 of the 5 skim-milk
series, to within 2–7 %. It must NOT be used for (a) ethyl pentanoate above 50 °C in water
or above 60 °C in skim milk, or (b) t-2-hexenal above 60 °C in skim milk — the three cells
the paper's own footnotes exclude. The isoamyl-acetate/skim-milk mid-range carries a
demonstrated ~1.8× error at 40 °C.** Everything else is good to ≈±10 %.

### 7.2 THE GRID — 5 compounds × 4 media × 6 temperatures. All [Z], dimensionless air/matrix.

Cells in **⚠️** are outside the validity window of their own ΔH fit (§6, footnotes ⁽¹⁾⁽²⁾).

**Isoamyl acetate**
| medium | 30 °C | 40 °C | 50 °C | 60 °C | 70 °C | 80 °C |
|---|---:|---:|---:|---:|---:|---:|
| water | 4.50e-3 [M] | 7.07e-3 | 1.08e-2 | 1.61e-2 | 2.35e-2 | 3.34e-2 |
| skim milk | 4.20e-3 [M] | 6.30e-3 | 9.22e-3 | 1.32e-2 | 1.85e-2 | 2.53e-2 |
| anhydrous milk fat | 5.70e-5 [M] | 1.04e-4 | 1.82e-4 | 3.08e-4 | 5.06e-4 | 8.08e-4 |
| full-fat cream | 1.40e-4 [M] | 2.36e-4 | 3.86e-4 | 6.12e-4 | 9.46e-4 | 1.42e-3 |
| **skim-milk retention %** | **6.7** | ⚠️10.9 (paper: 20) | 14.8 | 18.2 | 21.3 | **24.2** (paper: 26) |

**Amyl acetate**
| medium | 30 °C | 40 °C | 50 °C | 60 °C | 70 °C | 80 °C |
|---|---:|---:|---:|---:|---:|---:|
| water | 3.00e-3 [M] | 5.06e-3 | 8.25e-3 | 1.31e-2 | 2.02e-2 | 3.04e-2 |
| skim milk | 2.50e-3 [M] | 3.97e-3 | 6.13e-3 | 9.21e-3 | 1.35e-2 | 1.94e-2 |
| anhydrous milk fat | 1.40e-5 [M] | 2.46e-5 | 4.18e-5 | 6.86e-5 | 1.10e-4 | 1.71e-4 |
| full-fat cream | 6.00e-5 [M] | 1.10e-4 | 1.94e-4 | 3.30e-4 | 5.45e-4 | 8.75e-4 |
| **skim-milk retention %** | 16.7 | 21.5 | 25.8 | 29.5 | **29.5→32.9** (paper: ~30 at 70) | 36.0 |

**Ethyl pentanoate** — ⚠️ **water valid to 50 °C only; skim milk valid to 60 °C only**
| medium | 30 °C | 40 °C | 50 °C | ⚠️60 °C | ⚠️70 °C | ⚠️80 °C |
|---|---:|---:|---:|---:|---:|---:|
| water | 2.80e-3 [M] | 4.67e-3 | 7.53e-3 | ⚠️1.18e-2 | ⚠️1.81e-2 | ⚠️2.69e-2 |
| skim milk | 2.10e-3 [M] | 3.47e-3 | 5.57e-3 | 8.67e-3 | ⚠️1.32e-2 | ⚠️1.95e-2 |
| anhydrous milk fat | 1.30e-5 [M] | 2.63e-5 | 5.09e-5 | 9.48e-5 | 1.70e-4 | 2.95e-4 |
| full-fat cream | 3.70e-5 [M] | 7.57e-5 | 1.48e-4 | 2.79e-4 | 5.05e-4 | 8.84e-4 |
| **skim-milk retention %** | **25.0** (paper: 25 ✅) | 25.6 | 26.1 | ⚠️26.6 | ⚠️27.1 | ⚠️27.5 (paper: 48 ✗) |

**Hexanal** — fully valid across the range
| medium | 30 °C | 40 °C | 50 °C | 60 °C | 70 °C | 80 °C |
|---|---:|---:|---:|---:|---:|---:|
| **water** | **2.50e-3** [M] | **4.17e-3** | **6.74e-3** | **1.06e-2** | **1.62e-2** | **2.42e-2** |
| **skim milk** | **1.80e-3** [M] | **2.84e-3** | **4.36e-3** | **6.51e-3** | **9.51e-3** | **1.36e-2** |
| **anhydrous milk fat** | **8.20e-5** [M] | **1.33e-4** | **2.08e-4** | **3.18e-4** | **4.73e-4** | **6.89e-4** |
| **full-fat cream** | **1.60e-4** [M] | **2.55e-4** | **3.95e-4** | **5.96e-4** | **8.78e-4** | **1.26e-3** |
| **skim-milk retention %** | 28.0 | 31.9 | 35.4 | 38.5 (*"significantly different above 60 °C"*) | 41.3 | **43.8** (paper: 42 ✅) |

**t-2-hexenal** — ⚠️ **skim milk valid to 60 °C only; the whole skim-milk series measured on a browning system (§4)**
| medium | 30 °C | 40 °C | 50 °C | 60 °C | ⚠️70 °C | ⚠️80 °C |
|---|---:|---:|---:|---:|---:|---:|
| water | 4.40e-4 [M] | 8.00e-4 | 1.40e-3 | 2.38e-3 | 3.90e-3 | 6.24e-3 |
| skim milk | 6.40e-5 [M] | 1.19e-4 | 2.14e-4 | 3.71e-4 | ⚠️6.23e-4 | ⚠️1.02e-3 |
| anhydrous milk fat | 1.00e-5 [M] | 1.78e-5 | 3.06e-5 | 5.08e-5 | 8.20e-5 | 1.29e-4 |
| full-fat cream | 1.80e-5 [M] | 3.63e-5 | 7.00e-5 | 1.30e-4 | 2.32e-4 | 4.02e-4 |
| **skim-milk retention %** | **85.5** | 85.1 | 84.7 | 84.4 | 84.0 | 83.7 — **flat, as the paper states** |

### 7.3 ⚠️ THE PAPER'S OWN MEASURED AIR–WATER K IS 6.2× BELOW THE HENRY'S CONSTANT IT PRINTS [Z]

The Henry's constants in Table I convert to a dimensionless air/water partition coefficient
as `K_aw = H/(R_gas·T)` with `R_gas = 8.206 × 10⁻⁵ atm·m³·mol⁻¹·K⁻¹` at 298.15 K. Back-
extrapolating Table III to 25 °C with Table IV's water ΔH:

| compound | **Table I Henry ⇒ K_aw(25 °C)** [C] | **Meynier measured ⇒ K_aw(25 °C)** [Z from M] | **ratio** |
|---|---:|---:|---:|
| isoamyl acetate | 2.399 × 10⁻² | 3.549 × 10⁻³ | **6.76×** |
| amyl acetate | 1.586 × 10⁻² | 2.281 × 10⁻³ | **6.95×** |
| ethyl pentanoate | 1.520 × 10⁻² | 2.141 × 10⁻³ | **7.10×** |
| **hexanal** | **8.706 × 10⁻³** | **1.911 × 10⁻³** | **4.56×** |
| **t-2-hexenal** | 1.999 × 10⁻³ | 3.214 × 10⁻⁴ | **6.22×** |
| | | **geometric mean** | **6.24× (range 4.56–7.10)** |

**All five compounds, both chemical classes, the same direction, within a 1.6× spread. That
is a systematic calibration offset, not chemistry.** The most likely cause is named in §2:
**the FID was calibrated by liquid injections of aroma in cyclohexane and then used to
quantify a gas-phase headspace aliquot of 75–400 µL** — a liquid-standard-to-gas-sample
transfer with no gas-phase standard anywhere in the method.

**⚠️ The paper's own cross-lab comparison points the same way and the authors' explanation
does not cover it** (p. 227): *"air-water partition coefficients of 1.8 × 10⁻² and
3.8 × 10⁻³ at 30 °C for hexanal and t-2-hexenal were previously measured by Hall and
Anderson [16] against 2.5 × 10⁻³ and 0.44 × 10⁻³ in this study. The discrepancy between the
results can be attributed to the use of ammonium sulphate in the latter study."*
**Hall & Anderson are 7.2× and 8.6× HIGHER.** Salting-out with ammonium sulphate does raise
K, but the offset is the **same sign and similar size as the Table I discrepancy, where
there is no salt at all** — so salting-out cannot be the whole story.

**The corpus-level consequence, and it is directly relevant to the repo's hexanal
over-prediction:**

| source | hexanal air/water K, ~25–30 °C | provenance |
|---|---:|---|
| Meynier Table I (EPI Suite Henry) | **8.7 × 10⁻³** @25 °C | estimator [C] |
| **Meynier measured** | **2.5 × 10⁻³** @30 °C (**1.9 × 10⁻³** @25 °C) | **static HS, liquid-calibrated FID [M]** |
| Hall & Anderson 1983, via Meynier | **1.8 × 10⁻²** @30 °C | + ammonium sulphate [C] |
| **spread across the three** | **9.5×** | — |

**⇒ The literature range on the single most important partition constant in the repo's
hexanal lane is 9.5× wide, and the repo's likely default (an EPI-Suite/Buttery-class Henry
value ≈8.7 × 10⁻³) sits at 4.6× ABOVE the only direct measurement in this paper.** If the
repo computes matrix headspace ppb by multiplying a liquid concentration by an estimator
Henry constant, that alone contributes up to **4.6× of over-prediction** before any binding
or matrix term is applied. **This is a concrete, quantified, previously unrecorded candidate
for part of the known hexanal over-prediction — and, being a partition term, it is the
mechanism `k2_matrix_and_thresholds.md` §D.4.5 said the binding model *should* be doing its
job on.** ⚠️ **Not a licence to swap the constant**: the 6.2× is systematic across five
compounds and is at least as likely to be a Meynier method bias as a defect in the
estimator. **The finding is that the uncertainty on air/water K is ±0.5 decades, not that
Meynier is right.**

---

## 8. THE EMULSION MODEL — eqs. (4)–(6), and where it fails. [F]/[M]

Model as printed (p. 232), citing **Buttery, Guadagni & Ling 1973** [C]:

```
(4)  K_a,em = 1 / [ (Φ_cp / K_a,cp) + (Φ_dp / K_a,dp) ]
(5)  K_a,ffc = 1 / [ (0.7 / K_a,w)  + (0.3 / K_a,amf) ]     — water as continuous phase
(6)  K_a,ffc = 1 / [ (0.7 / K_a,sm) + (0.3 / K_a,amf) ]     — skim milk as continuous phase
```
with **Φ_dp = 0.3** the measured cream fat fraction. *"Equations (4) to (6) did not take into
account the presence of an interface layer."*

### 8.1 Measured vs calculated cream — the paper's findings, verbatim and complete

| compound | direction of failure | where it starts | the paper's wording |
|---|---|---|---|
| **isoamyl acetate** | **measured < calculated** (more retained than predicted) | **above 60 °C** | *"exhibited experimental coefficients significantly lower than the calculated ones above 60 ... °C"* |
| **hexanal** | **measured < calculated** | **above 40 °C** | *"...and 40 °C respectively"* |
| **t-2-hexenal** | **measured < calculated** (vs eq. 6) | **from 70 °C** | *"significantly lower than the calculated ones (Eq. (6) from 70 °C)"* |
| **amyl acetate** | **measured > calculated** (LESS retained than predicted) | **from 40 °C**, and the gap **grows with temperature** | *"exhibited a significantly greater experimental coefficient than the calculated one from 40 °C, coinciding with a higher release of amyl acetate from full-fat cream than expected, i.e. the presence of an interface increased the volatility of amyl acetate."* |
| **ethyl pentanoate** | **no significant difference, 30–80 °C** | — | *"presumably because the interface layer did not significantly influence the volatility of this aroma."* |

**⇒ A two-phase volume-weighted partition model applied to a real 30 % o/w emulsion fails in
BOTH directions on chemically similar molecules — amyl acetate and isoamyl acetate are
structural isomers with log P 2.3 and 2.26, and the model over-predicts retention for one
and under-predicts it for the other.** The residual is assigned to the **interfacial protein
layer**, which the paper did not characterise. **This is a hard, measured bound on how well
any composition-only partition model can do in an emulsified matrix, and it is the
partition-side analogue of the Hong 2020 result that log P does not predict matrix threshold
shifts.**

Choice of continuous phase (p. 233): *"the nature of the continuous phase had no significant
effect on the calculated partition coefficient of the three esters"*; for hexanal eqs. (5)
and (6) agree **except at 80 °C**; **for t-2-hexenal the eq. (6) values are significantly
lower than eq. (5) from 50 °C to 80 °C** — *"for t-2-hexenal, the nature of the continuous
phase of the emulsion may be predicted to affect its partition behaviour in emulsion"*, and
*"the retention of t-2-hexenal by milk proteins was still visible in the presence of 30 %
lipids, indicating a strong interaction."*
Bottom line as printed: *"The simple model (Eq. (4)) ... fitted well with measured values for
temperatures below 40 °C. However, above this temperature, a discrepancy appeared."*

---

## 9. THE CITED BINDING AND RETENTION FACTS — all [C], useful, none measured here

Collected from the Discussion (pp. 231–232). These extend
`k2_matrix_and_thresholds.md` §(b) and §B.3 and are recorded with their original sources.

| quantity | value | system | cited source |
|---|---:|---|---|
| β-lactoglobulin affinity constant, **ethyl pentanoate** | **366** (units not given; M⁻¹ implied) | 30 g/L protein, **pH 3** | Reiners, Nicklaus & Guichard 2000 [25] |
| β-lactoglobulin affinity constant, **isoamyl acetate** | **543–627** | " | Guichard & Langourieux 2000 [14] |
| retention of **isoamyl acetate** by β-lg | **1 to 16 %** | " | Guichard & Langourieux 2000 |
| retention of **2-nonanone** by β-lg | **44 to 60 %** | " | Guichard & Langourieux 2000 |
| **heptanal vs 2-nonanone binding to whey protein** | **heptanal > 2-nonanone at pH 6.9; the reverse at pH 4.7** | whey proteins | Mills & Solms 1984 [23] |
| **heptanal binding is "only partially reversible"** | qualitative | whey proteins | Mills & Solms 1984 |
| **retention of hexanal by 50 g/L soy protein** | **37 to 44 %** | soy protein solution | Gremli 1974 [13] |
| **retention of t-2-hexenal by 50 g/L soy protein** | **68 to 75 %** | " | Gremli 1974 |
| ethyl hexanoate retention by casein | **38 % at 5 g/L → 61 % at 50 g/L** | caseinate solutions | Landy, Druaux & Voilley 1995 [19] |
| ester retention by dairy proteins (30 g/L) | **6 % (ethyl butanoate) to 40 % (ethyl hexanoate)** | whey + β-lg mixtures | Fabre, Aubry & Guichard 2002 [11] |
| **2-octenal binds covalently to the imidazole ring of histidine** in BSA | qualitative mechanism | bovine serum albumin | Alaiz & Giron 1994 [1] |
| histidine-containing peptides quench **t-2-hexenal more than hexanal** | qualitative, direction only | model | Zhou & Decker 1999a,b [32,33] |
| UHT causes **partial irreversible unfolding of β-lg, increasing surface hydrophobicity** | mechanism for the observed rise of retention with temperature | — | Swaisgood 1996 [28] |

**⚠️ Gremli 1974's soy numbers (hexanal 37–44 %, t-2-hexenal 68–75 % at 50 g/L) are a
second, independent report of the same alkenal/alkanal asymmetry Meynier measures in skim
milk (1.39× vs 6.88×).** Two matrices, two decades apart, same direction. **Convert to a
per-gram basis for comparison with `k2_matrix_and_thresholds.md` §(b):** at 50 g/L, 40 %
retention ⇒ apparent `K_eff = 1.33 × 10⁻² L/g` for hexanal and, at 71 %, `4.90 × 10⁻² L/g`
for t-2-hexenal **[Z]**. The hexanal figure sits between Damodaran's dialysis 1.47 × 10⁻³
and Barallat-Pérez's headspace 5.15 × 10⁻² — **consistent with K2 §B.3's ruling that the
method (headspace counts irreversible loss, dialysis does not) is a first-class field**, and
Gremli is a **vacuum-extraction/headspace** method, so it belongs on the headspace side.

---

## 10. NAMED LAUNDERING HAZARDS

| # | claim, as printed | reality | anchor |
|---:|---|---|---|
| M-1 | Abstract: **"a significant retention of t-2-hexenal was observed in skim milk (nearly 90 % whatever the temperature)"** | Results (p. 228) say **"approximately 85 %"** and Table III gives **85.5 %**. The abstract rounds 85.5 up to "nearly 90". | Abstract vs p. 228 / T3 |
| M-2 | Abstract: **"the retention of the other aromas varied from 6 % for isoamyl acetate to 40 % for hexanal in skim milk"** | The **6 %** is isoamyl acetate at **30 °C**; the **40 %** is hexanal at **~78–80 °C**. Presented as a range, they are **two different temperatures**. At a common temperature the range is 6.7–28.0 % (30 °C) or 24.2–43.8 % (80 °C). | Abstract vs §7.2 |
| M-3 | **Table I** — MW, BP, P_sat, S, log P, Henry's constant and γ | **Every value is estimated/cited from the EPI-Suite demo server; the γ column is a UNIFAC calculation. Nothing in Table I was measured in this paper.** | T1 footnotes, p. 225 |
| M-4 | **The t-2-hexenal / skim-milk partition coefficient** presented alongside the other 19 | Measured at **100–5 800 mg/L, up to 1 933× the water arm**, on a system that **visibly browned**. The paper itself says the constant "must be considered as *apparent*" — but the caveat appears once, in prose, not in the table. | T2 + p. 228 |
| M-5 | The **ΔH values for ethyl pentanoate (water, skim milk) and t-2-hexenal (skim milk)** presented in the same table as the other 17 | Fitted over **30–50 °C** and **30–60 °C** only, per the table's own superscript footnotes ⁽¹⁾⁽²⁾. Extrapolating them to 80 °C is wrong by up to **1.7×** and the paper says the curve "levels off". Easy to miss: the footnote markers are 6 pt digits inside the cells. | T4 footnotes, p. 230 |
| M-6 | **"The discrepancy [with Hall & Anderson, 7.2–8.6×] can be attributed to the use of ammonium sulphate in the latter study."** | Salting-out is real, but the paper's own measured K is also **4.6–7.1× below the salt-free Henry's constants printed in its own Table I**. A single-cause explanation does not cover a systematic offset that appears against a no-salt reference too. | p. 227 vs §7.3 |
| M-7 | **FID "calibrated by manual injections of aroma compounds in cyclohexane"**, then used on gas-phase headspace aliquots | A **liquid-standard → gas-sample** calibration transfer with no gas standard. The reported inter-operator CV (<5 %) measures reproducibility, **not** accuracy, and cannot detect this class of error. | p. 226 |

**Plus two internal contradictions recorded but not resolvable from the PDF:** (a) the
skim-milk retention of isoamyl acetate at 40 °C — the prose says **20 %**, the paper's own
Table III + Table IV reconstruct **10.9 %**; (b) ethyl pentanoate skim-milk retention at
80 °C — prose says **48 %**, reconstruction says **27.5 %**, and the ΔH footnote says the fit
does not reach 80 °C, so **the prose value must come from the figure and cannot be
reproduced from any printed number.**

---

## 11. CONSOLIDATED NEW-PARAMETER TABLE

**Common conditions:** static headspace GC-FID, 3.000 g sample in a 22.4 mL vial, 30 min
equilibration (**60 min for skim milk**), K = C_gas/C_liquid dimensionless, mean of 10
(5 concentrations × duplicate), INRA Nantes 2002.

| # | parameter | value | units | matrix / conditions | class | anchor |
|---:|---|---:|---|---|:--:|---|
| 1–5 | **K_air/water, 30 °C** — isoamyl acetate / amyl acetate / ethyl pentanoate / **hexanal** / **t-2-hexenal** | **4.5e-3 ± 0.6e-3 / 3.0e-3 ± 0.4e-3 / 2.8e-3 ± 0.1e-3 / 2.5e-3 ± 0.1e-3 / 4.4e-4 ± 0.7e-4** | dimensionless | ultrapure water, 30 °C | M | T3 p.226 |
| 6–10 | **K_air/skim-milk, 30 °C**, same order | **4.2e-3 ± 0.5e-3 / 2.5e-3 ± 0.2e-3 / 2.1e-3 ± 0.3e-3 / 1.8e-3 ± 0.2e-3 / 6.4e-5 ± 2.7e-5** | dimensionless | UHT skim milk, <0.3 % fat, ~34 g/L protein (cited), 30 °C | M | T3 |
| 11–15 | **K_air/anhydrous-milk-fat, 30 °C**, same order | **5.7e-5 ± 0.3e-5 / 1.4e-5 ± 0.2e-5 / 1.3e-5 ± 0.2e-5 / 8.2e-5 ± 0.6e-5 / 1.0e-5 ± 0.1e-5** | dimensionless | 100 % milk fat, mp 32–35 °C, 30 °C | M | T3 |
| 16–20 | **K_air/full-fat-cream, 30 °C**, same order | **1.4e-4 ± 0.2e-4 / 6.0e-5 ± 0.3e-5 / 3.7e-5 ± 0.4e-5 / 1.6e-4 ± 0.3e-4 / 1.8e-5 ± 0.5e-5** | dimensionless | 30 % milk fat o/w emulsion, 30 °C | M | T3 |
| 21–40 | **ΔH_affinity, 5 compounds × 4 media** | see §6 table; **32.0 to 56.5 ± 0.4–4.5** | kJ·mol⁻¹ | van 't Hoff slope of ln K vs 1/T, 30–80 °C ⚠️ **3 cells restricted to 30–50 / 30–60 °C** | F | T4 p.230 |
| 41–45 | ΔH_vaporisation, pure compounds | 44.0 / 41.9 / 46.0 / **45.3** / **44.8** | kJ·mol⁻¹ | pure liquid | C | T4 |
| 46 | **water/skim-milk headspace suppression, hexanal** | **1.39** | × | 30 °C, ~34 g/L dairy protein | Z | §5.1 |
| 47 | **water/skim-milk headspace suppression, t-2-hexenal** | **6.88** | × | " | Z | §5.1 |
| 48 | water/skim-milk suppression, 3 esters | **1.07 / 1.20 / 1.33** | × | " | Z | §5.1 |
| 49 | **water/amf suppression, hexanal** | **30.5** | × | 30 °C, 100 % milk fat | Z | §5.1 |
| 50 | water/amf suppression, esters | **78.9 / 214 / 215** | × | " | Z | §5.1 |
| 51 | **skim-milk retention of hexanal vs water** | **28.0 / 31.9 / 35.4 / 38.5 / 41.3 / 43.8** | % at 30/40/50/60/70/80 °C | reconstructed; 80 °C validated against printed 42 % | Z | §7.2 |
| 52 | **skim-milk retention of t-2-hexenal vs water** | **85.5 / 85.1 / 84.7 / 84.4 / 84.0 / 83.7** | % at 30→80 °C — **flat** | " ⚠️ 70–80 °C outside ΔH validity | Z | §7.2 |
| 53 | **K(T) grid** | 120 values | dimensionless | 5 compounds × 4 media × 6 temperatures | **Z (reconstructed, ±10 % where validated)** | §7.2 |
| 54 | **systematic offset: Table I Henry ⇒ K_aw vs measured K_aw** | **6.24 (range 4.56–7.10)** | × over-statement by the estimator | 25 °C, 5 compounds | Z | §7.3 |
| 55 | **hexanal air/water K, literature spread** | **1.9e-3 (Meynier @25 °C) to 1.8e-2 (Hall & Anderson @30 °C) = 9.5×** | dimensionless | — | Z + C | §7.3 |
| 56 | emulsion two-phase model failure | **measured/calculated diverges in both directions above 40 °C; ethyl pentanoate alone is NS across 30–80 °C** | — | 30 % o/w cream | M | §8.1 |
| 57 | Gremli 1974 soy retention ⇒ per-gram | **hexanal 1.33e-2, t-2-hexenal 4.90e-2** | L/g | 50 g/L soy protein, vacuum extraction | Z from C | §9 |
| 58 | β-lg affinity constants (cited) | ethyl pentanoate **366**, isoamyl acetate **543–627** | M⁻¹ (implied) | 30 g/L, pH 3 | C | §9 |

---

## 12. PROPOSED FIT / HOLD-OUT ROLE — **DRAFT FOR ORCHESTRATOR**

> ⚠️ **Proposal only.** Meynier is a **new source**; a declaration amendment must be approved
> before any wave fits any row. This dossier does not edit the declaration.

| dataset | rows | **proposed role** | rationale |
|---|---:|---|---|
| **Table III, 20 K values at 30 °C** | 20 | **FIT-ELIGIBLE for RATIOS ONLY; HOLD-OUT for absolutes** | The 6.24× systematic offset (§7.3) makes the absolute scale untrustworthy but leaves the within-study matrix/water ratios intact — the offset cancels exactly in a ratio. **Propose: `K_matrix/K_water` ratios fittable; raw K values NOT.** |
| **Table IV, 20 ΔH_aff** | 20 | **FIT-ELIGIBLE (17 rows) / EXCLUDE (3 rows)** | Exclude ethyl pentanoate water + skim milk and t-2-hexenal skim milk (§6 footnotes ⁽¹⁾⁽²⁾). The 17 remaining are the paper's own regression outputs with SDs and are the cleanest temperature-dependence parameters in the wave. |
| **§7.2 reconstructed K(T) grid** | 120 | **DERIVED — usable as a PRIOR, never as data** | It is a reconstruction from rows 1–40, not an independent measurement. Using it as a fit target would double-count Tables III and IV. Mark `derived_from: Meynier_T3+T4`. |
| **The hexanal water column (2.5e-3 @30 °C → 2.42e-2 @80 °C)** | 6 | **HOLD-OUT — recommended as a direct probe of the over-prediction** | If the repo's hexanal partition term, evaluated at 30 °C in water, lands outside **1.9e-3 to 1.8e-2**, it is outside the entire measured literature. This is a cheap, decisive guard. |
| **§7.3 estimator-vs-measurement offset (6.24×)** | 1 | **HOLD-OUT / DIAGNOSTIC** | Not a parameter. It is the size of the uncertainty band that should sit on any air/water partition constant the repo uses (**±0.5 decades**). |
| **§8.1 emulsion model failure directions** | 5 | **HOLD-OUT — falsification test** | Any composition-only two-phase partition model the repo ships must be checked against the fact that this one fails in **both directions** on two structural isomers. |
| **Table I (7 columns × 5 compounds)** | 35 | **REJECT as measurements; INGEST only as `cited_from_EPI_Suite`** | Estimator output, not data (M-3). |
| **§9 cited binding constants** | 12 | **REJECT as parameters of this paper; RETAIN with original attribution** | All second-hand. The Gremli 1974 soy retention pair (row 57) is the one worth chasing to primary. |
| **t-2-hexenal / skim milk, all rows** | 8 | **QUARANTINE** | Measured on a browning system at up to 5 800 mg/L (M-4). Directionally sound, quantitatively unusable. |

---

## 13. RETRIEVALS THIS PAPER MAKES WORTH REQUESTING

1. **Hall, G. & Anderson, J. (1983)**, *Lebensm. Wiss. Technol.* **16**, 362–366 —
   *"Volatile fat oxidation products. II. Influence of temperature on volatility of
   saturated, mono- and di-unsaturated aldehydes in liquid media."* **The only other direct
   measurement of hexanal and t-2-hexenal air/water K vs temperature the corpus knows of,
   and it disagrees with Meynier by 7.2–8.6×.** Settling that disagreement is the highest-
   value retrieval for the hexanal lane. Also covers **di-unsaturated** aldehydes, i.e. it
   probably contains a t,t-2,4-decadienal partition value — the compound that is the extreme
   outlier in all three media of `k2_matrix_and_thresholds.md` §A.8 and for which the corpus
   has **no** partition constant at all.
2. **Gremli, H. A. (1974)**, *J. Amer. Oil Chem. Soc.* **51**, 95A–97A — *"Interaction of
   flavor compounds with soy protein."* The primary source of the hexanal 37–44 % /
   t-2-hexenal 68–75 % retention at 50 g/L soy protein (§9). **A soy-protein aldehyde
   retention pair is directly on the repo's plant-protein path** and would give a second
   headspace-method aldehyde constant to sit beside Barallat-Pérez's lupin values.
3. **Buttery, R. G., Guadagni, D. G. & Ling, L. C. (1973)**, *J. Agric. Food Chem.* **21**,
   198–201 — the source of the eq. (4) emulsion model **and** of the classical aqueous
   aldehyde partition/threshold values. Would let the repo's emulsion term and its aqueous
   Henry values be traced to one primary source.
4. **Landy, P., Druaux, C. & Voilley, A. (1995)**, *Food Chem.* **54**, 387–392 — casein
   retention **as a function of protein concentration (5 → 50 g/L)**. The corpus has almost
   no concentration-resolved retention data; this is a dose–response on protein loading.
