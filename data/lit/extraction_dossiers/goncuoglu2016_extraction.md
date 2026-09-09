# Göncüoğlu Taş & Gökmen 2016 — EXTRACTION (whole hazelnuts roasted at 150/160/170 °C for 15-120 min; 26-step multiresponse network on eleven measured responses; 26 rate constants at three temperatures with 95 % HPD, reference temperature T_b = 160 °C)

### A SECOND LABORATORY'S MULTIRESPONSE FIT IN A REAL FOOD: Hacettepe's hazelnut network prints all 26 constants at 150, 160 and 170 °C in a table that came through the text layer clean — and six of its steps are the same transformations the trunk carries, three of which (3-DG → 3,4-DG, 3,4-DG → HMF, glyoxal sink) agree with the shipped constants inside a factor of two while three others (1-DG → dimethylglyoxal 466x fast, HMF sink 23 000x fast, glucosone absent altogether) contradict them outright.

**Source on disk:** `data/articles/Goncouglu2016.pdf` (53 numbered pages in an Elsevier accepted-manuscript
PDF with line numbers, PII S0308-8146(16)32008-8).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/Goncouglu2016.txt`, 2160 lines). **Table 1 — the whole rate-constant table —
came through clean and is re-typed in full below.** Appendix A (the differential equations) is
garbled: every rate-constant subscript is lost in the text layer and the equations are NOT used
here; the network below is taken from Table 1's own step list, which names reactant and product for
all 26 steps, and from the Figure 2b line-art (which survives as ASCII and corroborates the step
list). **Supplementary Tables S1-S4 and Figures S1-S10 are not in the PDF.** In particular
**Table S4, the only place the authors' activation energies and reparametrised k_b values appear,
is absent from disk**, so *no activation energy from these authors exists in this repository* —
only the three temperature columns.

**Repo status before this dossier.** The same PDF already has two dossiers, written in earlier
waves: `goncuoglutas2016_extraction.md` (2026-09-07) and `goncuoglutas2017_extraction.md`
(2026-08-29, whose §0 records that the filename year is right and its own title year is wrong —
there is one paper, not two). Both are read-only extractions of this same accepted manuscript. This
dossier is the third read and the first written to the 2026-09-09 format, and its purpose is
different from theirs: to set Table 1 side by side with the constants the engine actually ships.
Nothing here contradicts either sibling; where a number is also in them, it is the same number.

## 0. Identity

| field | value |
|---|---|
| Title | "Maillard Reaction and Caramelization during Hazelnut Roasting: A multiresponse kinetic study" |
| Authors | Neslihan Göncüoğlu Taş, Vural Gökmen (corresponding, vgokmen@hacettepe.edu.tr) — Food Quality and Safety (FoQuS) Research Group, Department of Food Engineering, Hacettepe University, Beytepe Campus, 06800 Ankara, Turkey |
| Venue | Food Chemistry, accepted manuscript, reference FOCH 20292. Received 30 August 2016, revised 27 November 2016, accepted 30 November 2016 |
| DOI as printed | `http://dx.doi.org/10.1016/j.foodchem.2016.11.159` |
| Funding | TÜBİTAK IntenC (German-Turkish) project 113O178 |
| Naming | AP = Amadori product; HP = Heyns product; FFC = fructofuranosyl cation; 1,2-ED = 1,2-enediol; 1-DG / 3-DG / 3,4-DG = 1-deoxyglucosone / 3-deoxyglucosone / 3,4-dideoxyglucosone; GO = glyoxal; MGO = methylglyoxal; **DMG = dimethylglyoxal, i.e. 2,3-butanedione, i.e. diacetyl** (confirmed by the standards list — quinoxaline, 2-methylquinoxaline, 2,3-dimethylquinoxaline — and by the SIM ion 159); AA = **total** amino acids, one lumped pool; P1-P9 = unidentified product pools |
| Lineage | multiresponse machinery from van Boekel 1996 and Martins & van Boekel 2005 (both cited by name); the α-dicarbonyl derivatisation from Degen, Hellwig & Henle 2012; the HMF method from Kocadağlı et al. 2012; the sibling network from Kocadağlı & Gökmen 2016 (glucose/wheat flour) |
| Companions on disk | `goncuoglutas2016_extraction.md` and `goncuoglutas2017_extraction.md` (the same PDF); `kocadagli2016jafc_extraction.md` (the amine-free glucose glass that supplies the trunk's furanic and dicarbonyl constants) and `kocadagli2016foodchem_extraction.md` (glucose/wheat flour, same lab, same year, same dicarbonyl set); `martins2005_extraction.md` (the trunk itself) |

## 1. Why it matters

**The trunk is one laboratory's fit.** `src/kinetic_core/parameters.py`'s `MARTINS_M4` block is ten
constants from Martins & van Boekel 2005 — Wageningen, glucose 200 mmol/L + glycine 200 mmol/L in
0.1 mol/L phosphate at initial pH 6.8, aqueous, 80-120 °C — and `results/validation/kinetic_core_b1_fit_report.json`
is the fit built on it. The furanic channel (wave B7, `parameters_furanic.py`) and the dicarbonyl
trio (wave B13, `parameters_dicarbonyl.py`) both come from **Kocadağlı & Gökmen 2016 JAFC**, an
amine-free glucose *glass* at 160-200 °C with T_b = 180 °C.

**This paper is a real-food multiresponse fit from the Kocadağlı/Gökmen laboratory on eleven
responses at once.** Against Martins it is a second laboratory, a second matrix and a second
temperature window; against the B7/B13 constants it is the *same* laboratory in a *different*
matrix, which is a matrix test rather than a laboratory test, and this dossier says so at every row
rather than claiming more than it has.

Six steps in Table 1 are transformations the engine already carries as named reactions in
`src/kinetic_core/network.py`:

| this paper's step | Table 1 number | repository reaction (`network.py`) | registry key | wave |
|---|---|---|---|---|
| 3-DG → 3,4-DG | k18 | `r_tdg_ddg` | `k_tdg_ddg` | B7 |
| 3,4-DG → HMF | k19 | `r_ddg_hmf` | `k_ddg_hmf` | B7 |
| 1-DG → DMG (diacetyl) | k17 | `r_odg_da` | `k_odg_da` | B13 |
| GO → P6 (the glyoxal sink) | k25 | `r_go_sink` | `k_go_sink` | B13 |
| HMF → P7 (the HMF sink) | k26 | `r_hmf_self` | `k_hmf_self` | B7 |
| AP → 1-DG | k15 | `r_ama_odg` | `k_ama_odg` | B1 (Martins step 7) |
| AP → 3-DG | k14 | `r_ama_tdg` | `k_ama_tdg` | B1 (Martins step 4) |
| GLC + AA → AP | k5 | `r_schiff` + `r_amadori` (the composite) | `k_schiff` | B1 (Martins step 1) |

and three more bear directly on open questions the pre-registrations record:

1. **`kinetic_core_b13_prereg.md` §5 names, as the first item on the wishlist, "a glyoxal loss rate
   at two temperatures (a barrier for the sink)".** The trunk's `k_go_sink` carries
   Kocadağlı's rate with the **barrier fixed to zero by those authors**, so it runs at its 180 °C
   value at every temperature — the reason B13's own outcome section records the unexpected result
   that glucosone *accumulates above* glyoxal at 120 °C. **This paper prints that glyoxal loss rate
   at three temperatures: 18, 61 and 290 (×10⁻³ min⁻¹) at 150, 160 and 170 °C.** It is a 16-fold
   rise over 20 °C, so the sink emphatically does not have a zero barrier, and the wishlist item is
   answered in a real food (Flags 4 on why the answer is not yet a constant).
2. **`kinetic_core_b21_prereg.md` §1 and §6 turn on where glyoxal comes from.** The trunk makes it
   only through glucosone (`r_glc_g` then `r_g_go`), a route whose glass constants are a millionth
   per minute at 100 °C, and B21 fixed that by feeding glucosone from the Amadori compound on milk
   constants — at the cost, recorded in §6's "unforeseen finding", of worsening Martins' own Amadori
   series by 24 % in half-sum-of-squares. **This paper tested both routes and kept neither of them.**
   Its model discrimination retains **GLC → GO directly (k6)**, and the Results state plainly:
   "**Glucosone was not present in hazelnuts roasted at selected time-temperature combinations.**"
   That is a third topology for the glyoxal supply, from a laboratory that measures glucosone
   routinely (it is one of the seven quinoxalines in their SIM list, ion 251) and did not find it.
3. **`k5a_hmf_synthesis.md` gap G2, carried into `k_hmf_self`'s flags, is that the corpus has no HMF
   sink between 50 and 150 °C** — the shipped constant is 8.97×10⁻⁷ min⁻¹ derived from 0.9 % loss in
   seven days at **5 °C**, with Ea = 0 by declaration, and `parameters_furanic.py` pre-registers the
   consequence that HMF must be over-predicted. **This paper prints an HMF sink at 150, 160 and
   170 °C: 12, 21 and 103 (×10⁻³ min⁻¹).** At 160 °C that is **23 000 times** the shipped value.
   The pre-registered over-prediction now has a measured size.

What this paper does **not** give the repository: any activation energy (Table S4 absent); any
concentration-time datum as a number (all of Figure 1 is figure-only); any pH or moisture value
(Table S1 absent); any colour or browning number (Table S3 absent); any individual amino acid
(Table S2 absent); any melanoidin measurement at all; any water-activity measurement of its own.

## 2. Methods as they matter to a model

- **Matrix.** Whole unshelled hazelnuts (*Corylus avellana* L.), Turkish variety **Tombul**, from a
  local manufacturer in Giresun. **Five grams** per run. This is a *real food*: >56 % oil,
  11.7-20.8 % protein, and the sugars and amino acids that react are a few per cent of the mass.
- **Water.** The paper measures moisture (AOAC 925.10, 2 g dried at 105 °C to constant weight) and
  reports it **only in Table S1, which is absent**. The text says the moisture content decreased
  during roasting and that raw hazelnuts are **2.49-5.25 %** moisture with **water activity 0.40 to
  0.55** — but both of those are quoted from the literature (Köksal 2006; Lopez 1995), **not
  measured here**. Treat a_w 0.40-0.55 as a matrix description, not as a measurement of these
  samples, and note that it falls as the nuts roast.
- **pH.** Measured (0.5 g ground nut in 10 mL water, shaken 10 min, centrifuged, MeterLab PHM210)
  and reported **only in the absent Table S1**. The text says it decreased during roasting. **No pH
  number is on disk.** Every constant below is therefore at an unrecorded, drifting pH.
- **Heating.** Memmert UNE 400 oven at **150, 160 and 170 °C** for **15, 30, 60, 90 and 120 min**;
  nuts moved to a −18 °C freezer immediately on removal; skins removed by hand at room temperature;
  ground and kept at −18 °C. **No thermocouple.** The authors argue the heat-up is negligible
  because the specific heat of hazelnut (1994 J/kg·K, from Demir 2002) is about that of oil and half
  that of water, and cite Demir's own measurement that hazelnuts take **4-6 min** to reach oven
  temperature — then state that in their hands it "could take less time ... as there were no
  thermocouples incorporated". **Heat and mass transfer coefficients were deliberately not put into
  the model**, and the fit is isothermal from t = 0. The earliest sample is at 15 min, so the
  heat-up is a smaller fraction of the shortest run than in a typical aqueous pot, but it is not
  zero (Flags 3).
- **Extraction.** 0.50 g ground nut, triple water extraction (5 + 2.5 + 2.5 mL = 10 mL), vortex
  3 min and 5000 g for 5 min at each step, supernatants pooled and re-centrifuged 3 min.
- **Sugars.** Carrez I/II clarification, Oasis HLB cartridge, Agilent 1100 HPLC-RI, Shodex RSpak
  KC-811 (300 × 8 mm, 7 µm) at 40 °C, 0.1 % H₃PO₄ in water at 1.2 mL/min, 10 µL injection.
  External calibration **0.1-1 g/100 mL** for sucrose, glucose and fructose.
- **HMF.** Method of Kocadağlı et al. 2012 after 0.45 µm nylon filtration; external calibration
  **1-10 mg/L**.
- **α-Dicarbonyls.** o-Phenylenediamine derivatisation after Degen, Hellwig & Henle 2012: 200 µL
  extract + 800 µL acetonitrile:water (5:3), 5000 g for 5 min; then 0.5 mL supernatant + 150 µL
  0.5 M sodium phosphate buffer **pH 7** + 150 µL 0.2 % o-PDA in 10 mM DETAPAC; filtered at once
  and held **2 h at room temperature in the dark**. HPLC-ESI-MS, Agilent 1200 + 6130 single
  quadrupole, SIM on [M+H]⁺ of the quinoxalines: **251 glucosone, 235 both 1- and 3-deoxyglucosone,
  217 3,4-dideoxyglucosone, 159 dimethylglyoxal, 145 methylglyoxal, 131 glyoxal.**
  Calibration: quinoxaline / 2-methylquinoxaline / 2,3-dimethylquinoxaline at **0.1-2 mg/L**;
  3-deoxyglucosone derivatised as a standard at **0.1-5 mg/L**. **1-DG, 3,4-DG and glucosone have no
  standard of their own** — they are quantified against 3-DG's curve (they share ion 235 with 3-DG
  or are semi-quantitated), which is the same semi-quantitation caveat `k_tdg_ddg` and `k_ddg_hmf`
  already carry from the sibling paper.
- **Free amino acids and protein-bound lysine.** Free: extract + equal volume acetonitrile, 5000 g
  3 min, 0.45 µm. Bound lysine: **50 mg ground nut + 10 mL 8 N HCl, headspace nitrogen-filled,
  110 °C for 23 h**, filtered, 100 µL dried under nitrogen and taken up in 1 mL water:acetonitrile
  (50:50). Both by Waters Acquity UPLC-ESI-MS/MS (positive mode, triple quadrupole), Atlantis HILIC
  150 × 2.1 mm 3 µm, gradient of 0.1 % formic acid in water (A) / in acetonitrile (B), 0.4 mL/min,
  15 % A → 40 % A in 4 min, held 3 min, back to 15 % in 1 min; column 40 °C, autosampler 10 °C;
  capillary 3.5 kV, cone 20 V, extractor 3 V, source 120 °C, desolvation 350 °C, N₂ at 900 L/h.
  External calibration **0.05-2.0 mg/L** for every amino acid.
- **Colour.** Computer-vision image analysis (Mogol & Gökmen 2014), L*, a*, b* — **in the absent
  Table S3**. No browning number is on disk.
- **Initial charge, as printed in the Results.** Sucrose **5.5 ± 0.1 g/100 g dry weight**; fructose
  **0.4 ± 0.02 g/100 g dw**; glucose **0.2 ± 0.06 g/100 g dw**; total free amino acids
  **2112 ± 49 mg/kg dw**; protein-bound lysine **5401 ± 50 mg/kg dw**; "the concentration of total
  amino acids was **7513 ± 87 mg/kg dw**" — which is the pool the model calls AA, and it is
  the free pool plus the bound lysine (2112 + 5401 = 7513, exactly; **mine**). **In mmol/kg dry
  weight (mine, using M = 342.30 for sucrose, 180.16 for the hexoses): sucrose 161 mmol/kg,
  fructose 22.2 mmol/kg, glucose 11.1 mmol/kg.** The figure axes are in µmol/kg and the printed
  sucrose axis reaches 200 000 µmol/kg, i.e. 200 mmol/kg, consistent with that conversion. A
  **mmol/L is not derivable** — the matrix is a 3-5 % moisture solid and no volume, density or
  water content is on disk (Flags 1).
- **Fitting.** Athena Visual Studio **version 14.2**; ordinary differential equations, one per
  species (Appendix A), numerically integrated; **non-linear regression with the determinant
  criterion** (van Boekel 1996); goodness of fit and **95 % highest posterior density intervals**
  used to judge the models. Model discrimination on three questions: (i) whether the glucose-fructose
  isomerisation passes through 1,2-enediol, (ii) whether sucrose degradation and HMF formation run
  through the fructofuranosyl cation, (iii) whether the α-dicarbonyls react with amino acids.
- **Reference temperature.** The reparametrised Arrhenius equation is used with
  **T_b = 160 °C = 433.15 K exactly**, and "the data at 150, 160, and 170 °C were fitted to
  experimental data all at once". So the 160 °C column of Table 1 **is** the estimated k_b. This is
  the same T_b as Knol 2005's T_av and is **not** the repository's T_REF_K = 373.15 K (100 °C), nor
  the 180 °C of the Kocadağlı glass the B7/B13 constants come from. See section 4 for what that
  costs.
- **Fitted temperature window: 150-170 °C, twenty degrees wide.** The authors say so themselves:
  "these findings represented a temperature range of 150-170 °C and it might be better to study a
  wider temperature range, by including lower temperatures (100-120 °C), to explain the Arrhenius
  behavior of the reactions."

## 3. Tables re-typed

### Table 1 (p. 52 of the manuscript, lines 854-860). "Reaction rate constants with 95 % highest posterior density (HPD) intervals at different temperatures according to the proposed kinetic model in Figure 2b for Maillard reaction and caramelization during roasting of hazelnuts."

Units are as printed in the table's own "rate constant" column: `min⁻¹ × 10³` for every first-order
step and `kg × µmol⁻¹ × min⁻¹ × 10³` for the three second-order steps (5, 9, 10). The `× 10³`
means the printed number is 10³ × k, i.e. **k = (printed number) × 10⁻³**.

| # | elementary reaction step | rate constant unit | k @ 150 °C | HPD | k @ 160 °C | HPD | k @ 170 °C | HPD |
|---:|---|---|---:|---|---:|---|---:|---|
| 1 | SUC → GLC + FFC | min⁻¹ × 10³ | 6.9 | ±0.8 | 15 | ±3.1 | 22 | ±2.0 |
| 2 | GLC → 1,2-ED | min⁻¹ × 10³ | 141 | ±27.0 | 473 | ±131 | 698 | ±178 |
| 3 | 1,2-ED → GLC | min⁻¹ × 10³ | 0 | ±0 | 8.5 | ±3.6 | 28 | ±8.3 |
| 4 | GLC → 3-DG | min⁻¹ × 10³ | 0.03 | ±0.02 | 0 | ±0 | 0 | ±0 |
| 5 | GLC + AA → AP | kg × µmol⁻¹ × min⁻¹ × 10³ | 0.0009 | ±0.0007 | 0.003 | ±0.001 | 0.009 | ±0.001 |
| 6 | GLC → GO | min⁻¹ × 10³ | 0.6 | ±0.2 | 2.5 | ±0.9 | 9.2 | ±1.0 |
| 7 | 1,2-ED → FRU | min⁻¹ × 10³ | 1.3 | ±0.9 | 1.8 | ±0.7 | 4.2 | ±1.9 |
| 8 | FRU → 1,2-ED | min⁻¹ × 10³ | 0 | ±0 | 0 | ±0 | 41 | ±14 |
| 9 | FRU + AA → HP | kg × µmol⁻¹ × min⁻¹ × 10³ | 0.00023 | ±0.00007 | 0.00062 | ±0.00015 | 0 | ±0 |
| 10 | FFC + AA → HP | kg × µmol⁻¹ × min⁻¹ × 10³ | 0.00094 | ±0.00028 | 0.00027 | ±0.00030 | 0.00004 | ±0.00002 |
| 11 | FFC → HMF | min⁻¹ × 10³ | 0.58 | ±0.13 | 0.57 | ±0.19 | 2.02 | ±1.32 |
| 12 | HP → 1-DG | min⁻¹ × 10³ | 0.23 | ±0.18 | 0.85 | ±0.58 | 267 | ±177 |
| 13 | HP → 3-DG | min⁻¹ × 10³ | 0.009 | ±0.004 | 0.022 | ±0.030 | 12 | **ind** ᵇ |
| 14 | AP → 3-DG | min⁻¹ × 10³ | 0 | ±0 | 0.62 | ±0.03 | 0 | ±0 |
| 15 | AP → 1-DG | min⁻¹ × 10³ | 3.47 | ±3.12 | 3.51 | ±1.33 | 0.56 | ±0.58 |
| 16 | 1-DG → MGO | min⁻¹ × 10³ | 7012 | ±5510 | 33016 | **ind** ᵇ | 47920 | ±33240 |
| 17 | 1-DG → DMG | min⁻¹ × 10³ | 371 | ±241 | 895 | ±581 | 1073 | ±618 |
| 18 | 3-DG → 3,4-DG | min⁻¹ × 10³ | 4.27 | ±0.57 | 29.4 | ±26.1 | 88.1 | ±22.7 |
| 19 | 3,4-DG → HMF | min⁻¹ × 10³ | 0 | ±0 | 134 | ±127 | 390 | ±111 |
| 20 | HP → P1 | min⁻¹ × 10³ | 11 | ±1.4 | 4.7 | ±5.6 | 59 | **ind** ᵇ |
| 21 | AP → P2 | min⁻¹ × 10³ | 140 | ±122 | 21.2 | ±34.4 | 5.24 | ±1.75 |
| 22 | 1-DG → P3 | min⁻¹ × 10³ | 122 | **ind** ᵇ | 0 | **ind** ᵇ | 404 | **ind** ᵇ |
| 23 | MGO → P4 | min⁻¹ × 10³ | 126 | ±106 | 737 | ±113 | 918 | ±639 |
| 24 | DMG → P5 | min⁻¹ × 10³ | 54 | ±41 | 130 | ±90 | 106 | ±65 |
| 25 | GO → P6 | min⁻¹ × 10³ | 18 | ±8.3 | 61 | ±25 | 290 | **ind** ᵇ |
| 26 | HMF → P7 | min⁻¹ × 10³ | 12 | ±3.7 | 21 | ±11 | 103 | ±63.7 |

Footnote ᵃ, verbatim: "SUC: sucrose; GLC: glucose; FRU: fructose; FFC: fructofuranosyl cation;
1,2-ED: 1,2-enediol; AP: Amadori product; HP: Heyns product; 1-DG: 1-deoxyglucosone; 3-DG:
3-deoxyglucosone; 3,4-DG: 3,4-dideoxyglucosone; GO: glyoxal; MGO: methylglyoxal; DMG:
dimethylglyoxal; HMF: 5-hydroxymethylfurfural; AA: total amino acids; P: products."

Footnote ᵇ, verbatim: "ind: indeterminate, which means a large uncertainty in the estimated
parameter within 95 % confidence interval."

**Seven cells are exactly 0 ± 0** (steps 3 @150, 4 @160 and @170, 8 @150 and @160, 9 @170, 14 @150
and @170, 19 @150, 22 @160). Those are estimates that went to the boundary, not measurements of
zero rate, and none of them may be read as "this reaction does not occur at this temperature" — the
paper's own discussion treats step 14 at 160 °C (0.62) as evidence that AP → 3-DG "was estimated at
only one roasting temperature". **Five cells are marked indeterminate**, and the text names exactly
those five: "reaction rate constants of some steps (k13, k16, k20, k22, k25) showed a higher
uncertainty and could not be estimated in this interval."

### THE BARRIER TABLE IS NOT ON DISK

Section 3.2 says: "The temperature dependence of elementary reactions during hazelnut roasting was
determined by the activation energies (Ea) and reaction rate constants (k_b) at reference
temperature of 160 °C (**Table S4 in Supplementary material**)." **Table S4 is not in the PDF.**
The only quantitative statement about the barriers that survives in the body text is:

> "The activation energies of elementary reaction steps were found to range between **0-1174
> kJ/mol** with **six zero** and a few relatively high values."

and the authors' own verdict, repeated in the abstract, the discussion and the conclusions:

> "The temperature dependence of the reactions was complicated and could not be explained by the
> Arrhenius equation."

Per the house rule that a constant fitted at a stated reference temperature must have that
temperature recorded: **the reference temperature is 160 °C and the barriers themselves are
missing.** No Ea from this paper may be entered anywhere. The 0-1174 kJ/mol range is a range, not a
value, and the "six zero" barriers are the authors' fixed-to-zero device — the same device
`k_ddg_hmf` already carries from the sibling glass paper.

### Every number printed in the running text

| quantity | value | where (manuscript line) |
|---|---|---|
| initial sucrose | 5.5 ± 0.1 g/100 g dw | Results 3.1, line 270 |
| initial fructose | 0.4 ± 0.02 g/100 g dw | line 274 |
| initial glucose | 0.2 ± 0.06 g/100 g dw | line 275 |
| sucrose degradation rate (= Table 1 step 1) | 0.0069, 0.015, 0.022 min⁻¹ at 150 / 160 / 170 °C | line 272 |
| sucrose loss at 120 min | 60 %, 70 %, 90 % at 150 / 160 / 170 °C | line 273 |
| fructose and glucose after 30 min @150 °C | 0.26 ± 0.02 and 0.08 ± 0.01 g/100 g dw | line 276 |
| fructose and glucose after 15 min @160 °C | 0.21 ± 0.01 and 0.07 ± 0.01 g/100 g dw | line 277 |
| fructose and glucose after 15 min @170 °C | 0.18 ± 0.01 and 0.06 ± 0.02 g/100 g dw | line 278 |
| total free amino acids, raw | 2112 ± 49 mg/kg dw | line 281 |
| protein-bound lysine, raw | 5401 ± 50 mg/kg dw | line 281 |
| total amino acids (the model's AA pool) | 7513 ± 87 mg/kg dw | line 284 |
| total amino acid loss at 120 min | 68 %, 81 %, 85 % at 150 / 160 / 170 °C | line 284 |
| **3-deoxyglucosone at 120 min** | **6.7 ± 0.1 mg/kg dw (150 °C); 6.1 ± 0.1 mg/kg dw (160 °C)** | line 304 |
| **3-deoxyglucosone maximum at 170 °C** | **5.4 ± 0.1 mg/kg dw at 60 min**, then falls | line 307 |
| 3,4-dideoxyglucosone | "approximately 5 times lower than those of 3-deoxyglucosone" | line 311 |
| **1-deoxyglucosone maxima** | **0.22 ± 0.01, 0.31 ± 0.03, 0.27 ± 0.01 mg/kg dw at 150 / 160 / 170 °C** | line 314 |
| **glyoxal in RAW hazelnut** | **1.7 ± 0.6 mg/kg dw**; rises "up to 4 times" after 15 min and then does not change at any roasting temperature | line 317 |
| **methylglyoxal maximum** | **6.6 ± 0.5 mg/kg dw at 160 °C, 90 min** | line 319 |
| dimethylglyoxal | forms at all temperatures; no significant difference (p < 0.05) between 150 and 160 °C at 120 min | line 320 |
| **glucosone** | **"Glucosone was not present in hazelnuts roasted at selected time-temperature combinations."** | line 266 |
| **mannose** | looked for, "could not be detected" | line 428 |
| **HMF at 120 min** | **104 ± 0.5, 238 ± 1.9, 278 ± 0.7 mg/kg dw at 150 / 160 / 170 °C** | line 332 |
| EFSA dietary HMF intake estimate | 1.6 mg/person/day | line 336 |
| mass balance, moles recovered @150 °C | 90 % at 15 min; 39 % at 120 min | line 348 |
| mass balance @160 °C | 71 % at 15 min; 29 % at 120 min | line 350 |
| mass balance @170 °C | 62 % at 15 min; 15 % at 120 min | line 350 |
| hazelnut specific heat (cited) | 1994 J/kg·K (Demir 2002) | line 367 |
| time to reach oven temperature (cited) | 4-6 min (Demir 2002) | line 370 |
| raw hazelnut moisture (cited) | 2.49-5.25 % (Köksal 2006) | line 45 |
| raw hazelnut water activity (cited) | 0.40-0.55 (Lopez 1995) | line 46 |
| methylglyoxal/glyoxal from heated olive oil (cited comparator) | 0.61 ± 0.03 mg/kg and ~0.5 mg/kg after 200 °C / 1 h (Fujioka & Shibamoto 2004) | line 325 |
| lipid share of the hazelnut's MGO and GO | "not expected to be more than 10 % of their total concentration" | line 561 |
| industrial hazelnut roasting practice (cited) | 100-160 °C, 10-60 min; typically 145 °C for 15 min | lines 55-58 |

**Every concentration-time course is FIGURE-ONLY.** Figure 1 gives eleven panels (sucrose, fructose,
glucose, HMF, 1-deoxyglucosone, methylglyoxal, dimethylglyoxal, glyoxal, total amino acids, and by
the caption also 3-DG and 3,4-DG) as µmol/kg hazelnut against roasting time, symbols for observed
and lines for predicted, colour-coded 150 / 160 / 170 °C. Figures 3, 4 and 5 are the three
model-discrimination overlays (fructose and glucose without 1,2-enediol; HMF without the
fructofuranosyl cation; amino acids with the dicarbonyl-amine reactions). Per house rule none of
them is typed as a number. What survives from the axes as a *bound*, not a datum: the sucrose panel
runs to 200 000 µmol/kg, fructose to 25 000, glucose to 12 000, HMF to 2000, 1-DG to 2.5, MGO to
120, DMG to 20, GO to 150, total amino acids to 60 000 µmol/kg.

### The three model-discrimination verdicts, as printed

1. **1,2-enolisation is required.** Omitting the 1,2-enediol makes both fructose and glucose
   "continuously decrease during roasting", which the data reject (Fig. 3). "1,2-enolization is one
   of the rate determining steps."
2. **HMF comes mainly from the fructofuranosyl cation, not from 3-deoxyglucosone.** Omitting the
   cation route leaves predicted HMF "far below the experimental values" (Fig. 4). Within the 3-DG
   route, **3,4-DG → HMF (k19) is fast and 3-DG → 3,4-DG (k18) is the rate-determining step**;
   the text puts k19/k18 at "almost 5 times".
3. **The α-dicarbonyls do NOT react appreciably with amino acids in this matrix.** Including the
   Strecker (binary dicarbonyl + amino acid) steps spoils the amino-acid fit (Fig. 5), so all of
   them were removed and only the unimolecular dicarbonyl sinks k22-k26 kept. Also removed:
   **3-DG and 3,4-DG sinks** ("their fits were better in that case"), **FRU → FFC (k27 = 0)**, and
   **MGO from glucose (k28) and from 3-DG (k29)** — only **1-DG → MGO** survived.

The printed ordering of the dicarbonyl sinks: **MGO (k23) > DMG (k24) > GO (k25) > HMF (k26)**, all
rising with temperature. At 160 °C the printed values are 737 > 130 > 61 > 21 (×10⁻³ min⁻¹), which
reproduces the stated ordering exactly.

### Arithmetic on the printed constants (all mine)

**1. Ratios inside the 160 °C column** (unit-free within a set of first-order constants).
k5/k9 = 0.003/0.00062 = **4.8**, the paper's "almost 5 times higher" for Amadori over Heyns
formation at 160 °C; at 150 °C it is 0.0009/0.00023 = **3.9**, and at 170 °C the Heyns route from
fructose is estimated at zero, so the ratio is undefined and the paper's claim that glucose
dominates "especially at higher roasting temperatures" rests on a boundary estimate.
k19/k18 at 160 °C = 134/29.4 = **4.6**, the paper's "almost 5 times"; at 170 °C it is 390/88.1 =
**4.4**. k16/k17 = 33016/895 = **37** at 160 °C (the paper says "around 20, 40 and 45-fold" at
150/160/170; my values are **18.9, 36.9, 44.7** — consistent).

**2. Three-point Arrhenius refits of Table 1 (mine, and they mostly fail).** With
T = 423.15 / 433.15 / 443.15 K and a least-squares line through ln k against 1/T:

| step | Ea from a 3-point refit (kJ/mol) | R² | comment |
|---|---:|---:|---|
| 1, SUC → GLC + FFC | 90.6 | 0.968 | plausible |
| 2, GLC → 1,2-ED | 125.2 | 0.926 | plausible, close to van Boekel's "most reactions ~120" |
| 5, GLC + AA → AP | 179.5 | 1.000 | high but well-determined by three points |
| 6, GLC → GO | 212.9 | 1.000 | very high |
| 16, 1-DG → MGO | 150.5 | 0.897 | the 160 °C point is indeterminate |
| 17, 1-DG → DMG | 83.2 | 0.882 | — |
| 18, 3-DG → 3,4-DG | 236.4 | 0.979 | **six times the sibling glass paper's 36.9 for the same step** |
| 23, MGO → P4 | 155.7 | 0.842 | — |
| 25, GO → P6 | 216.4 | 0.993 | the glyoxal sink; the glass paper fixed this barrier to **zero** |
| 26, HMF → P7 | 166.9 | 0.922 | the HMF sink; the repository carries **zero** |
| 12, HP → 1-DG | 547.4 | 0.875 | physically impossible as a barrier |
| 13, HP → 3-DG | 557.6 | 0.832 | physically impossible as a barrier |

**These refits are NOT admissible as barriers and are recorded only to show why.** They confirm the
authors' own verdict from their (absent) Table S4 — a 0-1174 kJ/mol spread — from the data on disk:
over a twenty-degree window with pools that are simultaneously being made, drained and lumped, the
three-point slope is not a barrier. Every one of them is classed `derived_assumption` in section 4
and none should enter a registry. What they *are* good for is a direction: **the sinks the
repository holds at Ea = 0 (glyoxal, HMF, 3,4-DG → HMF) all rise steeply with temperature here, so
Ea = 0 is wrong in sign-of-slope terms, not merely imprecise.**

**3. What "×10³" costs, checked against the text.** Step 1 prints 6.9 / 15 / 22 in the
`min⁻¹ × 10³` column, and the Results say "the degradation rates of sucrose ... were **0.0069,
0.015 and 0.022 min⁻¹**". The table's convention is therefore confirmed by the authors' own prose:
**k = printed × 10⁻³.** This is the one unit check the paper makes for us, and it fixes the
first-order column. It does not fix the second-order column (Flags 1).

**4. The mass balance is a warning, not a nuisance.** Only 39 %, 29 % and 15 % of the initial moles
of the eleven measured species remain at 120 min at 150, 160 and 170 °C. So **60-85 % of the carbon
this model tracks leaves through P1-P9, pools that are never measured.** Every sink constant
(k20-k26) is fitted to the *disappearance* of its parent, with nothing on the product side to
constrain it — exactly the caveat Knol 2005 states for its acrylamide sink. Read k22-k26 as loss
rates, never as named reactions.

**5. AA = 7513 is free + bound (mine).** 2112 (free) + 5401 (protein-bound lysine) = 7513 mg/kg dw,
the printed total to the digit. So the second-order constants k5, k9 and k10 are per unit of a pool
that is **72 % protein-bound lysine**, not free amino acid — a lumping the trunk's `k_schiff` (pure
glycine) and wave B20's glycation arm (which separates `LYSP` from `Gly` precisely because they are
not interchangeable) both refuse. Flags 2.

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Present and keyed: **`hmf`** and
**`2_3_butanedione`** (which is this paper's DMG / dimethylglyoxal / diacetyl). **Absent:**
glucose, fructose, sucrose, the fructofuranosyl cation, 1,2-enediol, the Amadori and Heyns products,
1-deoxyglucosone, 3-deoxyglucosone, 3,4-dideoxyglucosone, **glyoxal**, **methylglyoxal**,
**glucosone**, and every amino acid. The registry is a product/marker list; the network's species
names (`Glc`, `Fru`, `AMA`, `TDG`, `ODG`, `DDG`, `GO`, `MGO`, `G`, `DA`, `HMF`) are network-local
ids, not registry ids.

Every row below shares these conditions: **whole Tombul hazelnuts, 5 g, oven-roasted in air at 150 /
160 / 170 °C for 15-120 min; real food, >56 % oil, 2.49-5.25 % moisture and a_w 0.40-0.55 as
literature values for the raw nut and falling during roast (not measured here); pH measured but
reported only in the absent Table S1; no thermocouple and no heat-up correction; total amino acids
lumped into one pool AA that is 72 % protein-bound lysine; Athena Visual Studio 14.2, determinant
criterion, three temperatures fitted simultaneously; T_b = 160 °C = 433.15 K; barriers in the absent
Table S4.**

| step | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| SUC → GLC + FFC | k1 at 150/160/170 °C | 6.9 ±0.8 / **15 ±3.1** / 22 ±2.0 | ×10⁻³ min⁻¹ | as above | first order in sucrose | Table 1 step 1; confirmed in prose line 272 | measured_rate |
| GLC → 1,2-ED | k2 | 141 ±27.0 / **473 ±131** / 698 ±178 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 2 | measured_rate |
| 1,2-ED → GLC | k3 | 0 ±0 / **8.5 ±3.6** / 28 ±8.3 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 3 | measured_rate (150 °C cell on the boundary) |
| GLC → 3-DG | k4 | 0.03 ±0.02 / **0 ±0** / 0 ±0 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 4 | measured_rate (two cells on the boundary) |
| **GLC + AA → AP** | **k5** | 0.0009 ±0.0007 / **0.003 ±0.001** / 0.009 ±0.001 | **×10⁻³ kg·µmol⁻¹·min⁻¹** (= ×10⁻³ kg/(µmol·min); ×1 kg/(mmol·min) after the µmol→mmol factor, **mine**) | " | **second order**, glucose × total amino acids | Table 1 step 5 | measured_rate (unit basis is per kg of food, Flags 1) |
| GLC → GO | k6 | 0.6 ±0.2 / **2.5 ±0.9** / 9.2 ±1.0 | ×10⁻³ min⁻¹ | " | first order in glucose | Table 1 step 6 | measured_rate — **the surviving glyoxal route; no glucosone intermediate** |
| 1,2-ED → FRU | k7 | 1.3 ±0.9 / **1.8 ±0.7** / 4.2 ±1.9 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 7 | measured_rate |
| FRU → 1,2-ED | k8 | 0 ±0 / **0 ±0** / 41 ±14 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 8 | measured_rate (two cells on the boundary) |
| FRU + AA → HP | k9 | 0.00023 ±0.00007 / **0.00062 ±0.00015** / 0 ±0 | ×10⁻³ kg·µmol⁻¹·min⁻¹ | " | second order | Table 1 step 9 | measured_rate |
| FFC + AA → HP | k10 | 0.00094 ±0.00028 / **0.00027 ±0.00030** / 0.00004 ±0.00002 | ×10⁻³ kg·µmol⁻¹·min⁻¹ | " | second order | Table 1 step 10 | measured_rate (**falls with temperature**; the 160 °C HPD exceeds the estimate) |
| FFC → HMF | k11 | 0.58 ±0.13 / **0.57 ±0.19** / 2.02 ±1.32 | ×10⁻³ min⁻¹ | " | first order in an **unmeasured** cation pool | Table 1 step 11 | measured_rate — **ratio-only**: [FFC] is never measured, so only k10·[FFC] and k11·[FFC] are identified |
| HP → 1-DG | k12 | 0.23 ±0.18 / **0.85 ±0.58** / 267 ±177 | ×10⁻³ min⁻¹ | " | first order in an **unmeasured** Heyns pool | Table 1 step 12 | measured_rate — ratio-only; a 314x jump over 10 °C |
| HP → 3-DG | k13 | 0.009 ±0.004 / **0.022 ±0.030** / 12 (**ind**) | ×10⁻³ min⁻¹ | " | first order | Table 1 step 13 | measured_rate — **the authors name k13 as indeterminate**; the 160 °C HPD exceeds the estimate |
| **AP → 3-DG** | **k14** | 0 ±0 / **0.62 ±0.03** / 0 ±0 | ×10⁻³ min⁻¹ | " | first order in an **unmeasured** Amadori pool | Table 1 step 14 | measured_rate — one temperature only, the other two on the boundary |
| **AP → 1-DG** | **k15** | 3.47 ±3.12 / **3.51 ±1.33** / 0.56 ±0.58 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 15 | measured_rate — **falls at 170 °C**; the 170 °C HPD spans zero |
| 1-DG → MGO | k16 | 7012 ±5510 / **33016 (ind)** / 47920 ±33240 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 16 | measured_rate — **the largest constant in the paper (33 min⁻¹ at 160 °C) and the authors mark it indeterminate at T_b** |
| **1-DG → DMG (diacetyl)** | **k17** | 371 ±241 / **895 ±581** / 1073 ±618 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 17 | measured_rate |
| **3-DG → 3,4-DG** | **k18** | 4.27 ±0.57 / **29.4 ±26.1** / 88.1 ±22.7 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 18 | measured_rate (the 160 °C HPD is 89 % of the estimate) |
| **3,4-DG → HMF** | **k19** | 0 ±0 / **134 ±127** / 390 ±111 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 19 | measured_rate (150 °C on the boundary; 160 °C HPD 95 % of the estimate) |
| HP → P1 | k20 | 11 ±1.4 / **4.7 ±5.6** / 59 (**ind**) | ×10⁻³ min⁻¹ | " | first order | Table 1 step 20 | measured_rate — authors name k20 indeterminate; non-monotonic |
| AP → P2 | k21 | 140 ±122 / **21.2 ±34.4** / 5.24 ±1.75 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 21 | measured_rate — **falls monotonically with temperature** |
| 1-DG → P3 | k22 | 122 (**ind**) / **0 (ind)** / 404 (**ind**) | ×10⁻³ min⁻¹ | " | first order | Table 1 step 22 | measured_rate — **indeterminate at all three temperatures**; do not use |
| MGO → P4 | k23 | 126 ±106 / **737 ±113** / 918 ±639 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 23 | measured_rate |
| DMG → P5 (the diacetyl sink) | k24 | 54 ±41 / **130 ±90** / 106 ±65 | ×10⁻³ min⁻¹ | " | first order | Table 1 step 24 | measured_rate — non-monotonic |
| **GO → P6 (the glyoxal sink)** | **k25** | 18 ±8.3 / **61 ±25** / 290 (**ind**) | ×10⁻³ min⁻¹ | " | first order in glyoxal | Table 1 step 25 | measured_rate — authors name k25 indeterminate at 170 °C; the 150 and 160 °C cells are determinate |
| **HMF → P7 (the HMF sink)** | **k26** | 12 ±3.7 / **21 ±11** / 103 ±63.7 | ×10⁻³ min⁻¹ | " | first order in HMF | Table 1 step 26 | measured_rate |
| ratio, Amadori over Heyns formation | k5/k9 | 3.9 / **4.8** / undefined (k9 = 0) | — | " | — | derived from Table 1 (mine); the paper says "almost 5 times" | within_study_ratio |
| ratio, HMF step over its parent | k19/k18 | undefined (k19 = 0) / **4.6** / 4.4 | — | " | — | derived (mine); the paper says "almost 5 times" | within_study_ratio |
| ratio, MGO over DMG from 1-DG | k16/k17 | 18.9 / **36.9** / 44.7 | — | " | — | derived (mine); the paper says "around 20, 40 and 45-fold" | within_study_ratio |
| dicarbonyl sink ordering | k23 > k24 > k25 > k26 | 737 > 130 > 61 > 21 at 160 °C | ×10⁻³ min⁻¹ | " | — | Table 1 + text line 597 | within_study_ratio |
| 3-deoxyglucosone, roasted | 6.7 ± 0.1 (150 °C, 120 min); 6.1 ± 0.1 (160 °C, 120 min); 5.4 ± 0.1 (170 °C, 60 min, maximum) | mg/kg dw | " | — | Results line 304-307 | **level_only** |
| 1-deoxyglucosone, maxima | 0.22 ± 0.01 / 0.31 ± 0.03 / 0.27 ± 0.01 at 150 / 160 / 170 °C | mg/kg dw | " | — | Results line 314 | **level_only** |
| glyoxal, raw nut | 1.7 ± 0.6 | mg/kg dw | unroasted | — | Results line 317 | **level_only** — a background the model does not have |
| glyoxal, roasted | ~4x the raw value after 15 min, flat thereafter at every temperature (so ~6.8 mg/kg dw, **mine**, from "up to 4 times") | mg/kg dw | " | — | Results line 317 | **level_only** (the 4x is prose, not a table) |
| methylglyoxal, maximum | 6.6 ± 0.5 at 160 °C, 90 min | mg/kg dw | " | — | Results line 319 | **level_only** |
| **glucosone** | **not detected at any of the fifteen time-temperature combinations** | — | " | — | Results line 266 | **measured null** (the method quantifies it: SIM 251) |
| mannose | not detected | — | " | — | Discussion line 428 | measured null |
| HMF at 120 min | 104 ± 0.5 / 238 ± 1.9 / 278 ± 0.7 at 150 / 160 / 170 °C | mg/kg dw | " | — | Results line 332 | **level_only** |
| 3,4-dideoxyglucosone | "approximately 5 times lower" than 3-DG | — | " | — | Results line 311 | within_study_ratio (prose) |
| molar recovery of the eleven tracked species | 90/71/62 % at 15 min and 39/29/15 % at 120 min, at 150/160/170 °C | % of initial moles | " | — | Results line 347-351 | measured (a balance, not a rate) |
| initial charge | sucrose 161, fructose 22.2, glucose 11.1 | mmol/kg dry weight (**mine**, from the printed g/100 g dw) | raw nut | — | Results lines 270-275 | derived_assumption (arithmetic only) |
| total amino acid pool | 7513 ± 87 mg/kg dw = 2112 free + 5401 protein-bound lysine (**mine**, the sum is exact) | mg/kg dw | raw nut | — | Results lines 281-284 | measured, decomposition mine |
| three-point Arrhenius refits of Table 1 | see section 3 table (90.6 to 557.6 kJ/mol) | kJ/mol | " | — | derived from Table 1 (**mine**) | **derived_assumption — INADMISSIBLE as barriers**; recorded to show the temperature dependence is not Arrhenius |
| all eleven concentration-time courses at three temperatures | — | µmol/kg vs min | " | — | Figure 1 (and Figs. 3-5) | **figure_only** |
| pH and moisture against roasting time | — | — | " | — | Table S1, **absent from disk** | not on disk |
| individual free amino acids and bound lysine against time | — | mg/kg | " | — | Table S2, **absent** | not on disk |
| colour L*, a*, b* against time | — | — | " | — | Table S3, **absent** | not on disk |
| **activation energies and k_b at 160 °C** | — | kJ/mol | " | — | **Table S4, absent** | **not on disk — no barrier from this paper may be used** |

### Which of the trunk's constants can be set against this paper, and on what basis

Three things have to line up before a number here can be compared with a shipped one: **the unit
basis**, **the reference temperature**, and **the matrix**. They line up differently for each step,
so here is each one separately.

**The reference temperature is the easy part, once.** This paper's T_b is **160 °C = 433.15 K
exactly**, and the 160 °C column *is* k_b. The repository's own `T_REF_K` is **373.15 K (100 °C)**,
and the B7/B13 constants were re-referenced there from Kocadağlı's **180 °C**. So a comparison at
100 °C is impossible from this paper — it would need a barrier, and the barrier table is not on
disk. **Every comparison below is therefore made at 160 °C**, by transporting the *shipped* constant
up from 100 °C using the *shipped* barrier, and reading this paper's 160 °C column as printed. That
puts all the extrapolation risk on the repository's side, which is the honest place for it: the
shipped barrier is the thing under test.

**(a) `k_tdg_ddg`, 3-DG → 3,4-DG — comparable, and it AGREES.** Both first order, both min⁻¹, no
concentration basis needed. The shipped constant is Kocadağlı JAFC step 4: k_b = 30.5×10⁻³ min⁻¹ at
180 °C with Ea = 36.9 ± 6.3, which transports to **19.4×10⁻³ min⁻¹ at 160 °C (mine)**. This paper
prints **29.4 ± 26.1 ×10⁻³ min⁻¹ at 160 °C**. **A factor of 1.5**, and the hazelnut's HPD covers the
glass value comfortably. An amine-free glucose glass and a roasting nut agree on the
rate-determining step of the 3-DG limb within experimental scatter. This is the strongest
cross-matrix agreement in the paper. **The barriers do not agree**: my three-point refit here gives
236 kJ/mol against the glass's 36.9 — a sixfold difference that says the glass barrier, which
`parameters_furanic.py` already flags with `plateau_caveat_2.1x_over_40C`, is too flat.

**(b) `k_ddg_hmf`, 3,4-DG → HMF — comparable, and it AGREES.** Shipped: 119×10⁻³ min⁻¹ with
**Ea fixed to zero by Kocadağlı's own authors**, so 119×10⁻³ at every temperature. This paper:
**134 ± 127 ×10⁻³ min⁻¹ at 160 °C**. **A factor of 1.13** — the closest match on disk between the
engine and an independent matrix. But the agreement is at one temperature only: this paper's 0 /
134 / 390 across 150 / 160 / 170 °C is a steep rise where the shipped constant is flat, so the two
agree at 160 °C by construction of that flatness and would separate by 3x at 170 °C. **Both papers'
authors nevertheless describe this step as fast relative to its parent** — "almost 5 times higher"
here, 2-5x in the glass — and that ratio, not the absolute value, is the transportable claim.

**(c) `k_go_sink`, glyoxal → unassigned — comparable, and the RATE agrees while the BARRIER is
refuted.** Shipped: 32.6×10⁻³ min⁻¹, **Ea fixed to zero by Kocadağlı**, at every temperature. This
paper: **61 ± 25 ×10⁻³ min⁻¹ at 160 °C** — a factor of **1.87**, so the *rate* survives the move
from glass to nut. What does not survive is the zero barrier: **18 → 61 → 290 ×10⁻³ min⁻¹ over
150 → 170 °C is a 16-fold rise in twenty degrees.** `kinetic_core_b13_prereg.md` §5 asks for exactly
this ("a glyoxal loss rate at two temperatures (a barrier for the sink)") and §5's outcome section
records that the fixed-zero sink is what makes glucosone accumulate above glyoxal at 120 °C. **This
paper supplies the two temperatures.** It does not supply a usable barrier — the three-point refit
is 216 kJ/mol, which is not credible — but it establishes the *sign and rough size* of the missing
temperature dependence, which is enough to retire the claim that Ea = 0 is a neutral choice. The
170 °C cell is marked indeterminate by the authors; the 150 and 160 °C cells are not, and a
two-point slope on those two alone gives **190 kJ/mol (mine)**, still not credible, which is itself
the finding: over twenty degrees in a drying nut the slope is not a barrier. **What should be
requested is the same measurement over a wider window.**

**(d) `k_hmf_self`, the HMF sink — comparable, and the shipped value is WRONG BY FOUR DECADES.**
Shipped: **8.97×10⁻⁷ min⁻¹**, derived from Hamzalıoğlu & Gökmen 2018's model-free control (0.9 %
loss in 7 days at **5 °C**, pH 3.5), with **Ea = 0 by declaration** because there is only one
temperature. This paper: **21 ± 11 ×10⁻³ min⁻¹ at 160 °C**, i.e. **2.3×10⁴ times larger (mine)**.
`k5a_hmf_synthesis.md` gap G2 records that the 50-150 °C window is empty and
`parameters_furanic.py` pre-registers that HMF must therefore be over-predicted. **The
over-prediction now has a measured size in a real food**, and it is large: with the shipped sink,
HMF's half-life is 1.5 years at any temperature; with this paper's, it is **33 min at 160 °C
(mine)** and 6.7 min at 170 °C. This is the single most consequential row in the dossier for the
engine's HMF answers, and the sink is measured in the same food matrix whose HMF levels the paper
also prints (104 / 238 / 278 mg/kg dw at 120 min), so the sink and the level are internally
consistent within this one study.

**(e) `k_odg_da`, 1-DG → diacetyl — comparable, and it DISAGREES by two and a half decades.**
Shipped: Kocadağlı JAFC step 12, k_b = 12.2×10⁻³ min⁻¹ at 180 °C with Ea = 150.8 ± 8.8 (the
steepest barrier on the trunk), transporting to **1.92×10⁻³ min⁻¹ at 160 °C (mine)**. This paper's
k17, the same transformation (1-DG → dimethylglyoxal = 2,3-butanedione = diacetyl):
**895 ± 581 ×10⁻³ min⁻¹ at 160 °C**. **A factor of 466 (mine).** The two HPDs do not come within two
decades of each other. Either the steep glass barrier is badly wrong when extrapolated down from
180 °C, or the roasting nut has a route to diacetyl the amine-free glass does not (the paper's own
suggestion: Hollnagel & Kroh's amino-catalysed rearrangement, "an enhanced formation of
dimethylglyoxal in the presence of amino compounds"). **`k_da_sink` is a separate and equally sharp
conflict**: the trunk carries diacetyl's loss at **0 ± 0** (Kocadağlı step 17, diacetyl accumulates
in the glass) as an explicit prediction the data may reject, and this paper measures
**130 ± 90 ×10⁻³ min⁻¹ at 160 °C** (k24). **The prediction is rejected.**

**(f) The glucosone route, `r_glc_g` + `r_g_go` — this paper REFUSES the topology.** The shipped
route makes glyoxal only through glucosone; at 160 °C the shipped entry `k_glc_g` runs at
**1.47×10⁻⁵ min⁻¹ (mine)** and `k_g_go` at **0.233 min⁻¹ (mine, glass value)**. This paper's
network has **no glucosone node at all**, because "glucosone was not present in hazelnuts roasted at
selected time-temperature combinations" — from a laboratory that quantifies glucosone by its own
quinoxaline (SIM 251) as a matter of routine. Its glyoxal comes straight from glucose, **k6 =
2.5 ± 0.9 ×10⁻³ min⁻¹ at 160 °C**, which the shipped network has no slot for. This is a third
answer to the question `kinetic_core_b21_prereg.md` was written about: B13 had glucose → glucosone →
glyoxal on glass constants; B21 replaced the entry with Amadori → glucosone on Hamzalıoğlu's milk
constants and paid for it with a 24 % rise in the half sum of squares on Martins' own Amadori
series; **this paper says that in a low-moisture real food neither entry is needed and the sugar
goes to glyoxal directly.** It does not settle B21's conflict — a roasting nut at 150-170 °C is not
Martins' aqueous pot at 80-120 °C, and the constant is not transportable to water — but it removes
the assumption that the glucosone node is obligatory, and it is the first *measured null* on
glucosone in the corpus. **B21b, the pre-registered joint fit, should carry a direct
glucose → glyoxal edge as a third candidate topology.**

**(g) `k_ama_tdg` and `k_ama_odg`, the Amadori branch — SAME TRANSFORMATION, NOT COMPARABLE IN
MAGNITUDE.** Martins' step 4 (Amadori → 3-DG + Gly, k_ref 1.11×10⁻² min⁻¹, Ea 97.0 at 100 °C) and
step 7 (Amadori → 1-DG + Gly, 1.57×10⁻² min⁻¹, Ea 107.0) have exact counterparts here in k14
(AP → 3-DG, 0.62×10⁻³ min⁻¹ at 160 °C) and k15 (AP → 1-DG, 3.51×10⁻³). Both are first order in
min⁻¹, so the unit lines up and no concentration basis is needed. **But the Amadori pool is
UNMEASURED here** — the paper says so: "Although reactions leading to formation of Amadori/Heyns
product were involved in the proposed model, Amadori/Heyns product could not be measured
experimentally" — whereas Martins measured his. A first-order constant fitted against an unmeasured
pool is identified only up to that pool's scale, so **k14 and k15 are ratio-only** and may be
compared with Martins only as a ratio. That ratio is informative: **k15/k14 = 5.7 at 160 °C
(mine)**, i.e. the Amadori compound goes to 1-DG 5.7 times more often than to 3-DG, against
Martins' **k_ama_odg/k_ama_tdg = 1.57e-2/1.11e-2 = 1.41 at 100 °C (mine)**. Both put 1-DG ahead;
the hazelnut puts it four times further ahead. That is a real, transportable, cross-laboratory
agreement on **direction** and a disagreement on **degree**, and the degree is confounded by
temperature (160 vs 100 °C) and by matrix.

**(h) `k_schiff`, the condensation — NOT COMPARABLE without an assumption the paper does not
license.** Martins' step 1 is `1.6e-5 L/(mmol·min)` at 100 °C. This paper's k5 is
`0.003 ×10⁻³ kg·µmol⁻¹·min⁻¹` at 160 °C, i.e. **3×10⁻³ kg/(mmol·min) (mine)**. The unit is *per
kilogram of hazelnut*, not per litre of solution, and converting one into the other needs a volume —
which for a 3-5 % moisture nut is neither printed nor derivable (the reaction is not happening in
the free water; there is barely any). On top of that, Martins' amine is **glycine**, one compound at
a known 200 mmol/L, while this paper's is **AA, a lumped pool 72 % of which is protein-bound
lysine**. Two different quantities with the same dimension. **Do not put these two numbers side by
side.** What *is* comparable is the within-study statement both papers make about the same
structural question — that the aldose route dominates the ketose route — and here the paper prints
k5/k9 = 4.8 at 160 °C for glucose over fructose.

**(i) The isomerisation, `k_glc_fru` / `k_fru_glc` — STRUCTURALLY DIFFERENT, so no magnitude
transfers.** Martins writes Glc ⇄ Fru as two direct steps (1.64×10⁻³ and 9.15×10⁻³ min⁻¹ at
100 °C). This paper inserts an explicit **1,2-enediol** and splits each direction in two
(k2/k3 on the glucose side, k7/k8 on the fructose side), and its model discrimination shows the
enediol is required — without it "the concentrations of both fructose and glucose were not
estimated well and continuously decreased". The enediol pool is unmeasured, so all four constants
are ratio-only. **The structural claim is what transports: an explicit enediol is needed in a
low-moisture food and is not needed in Martins' pot.** The sibling paper's authors say why in
words this dossier records without endorsing: in the absence of amino acids the glucose-to-fructose
conversion is fast enough that the enediol need not be modelled (Kocadağlı & Gökmen 2016).

**(j) What cannot be transported at all.** No barrier (Table S4 absent). No pH (Table S1 absent).
No water activity of these samples. No browning or melanoidin measurement of any kind — so this
paper does **not** touch the trunk's browning hold-out, which remains the B1 fit's one out-of-sample
success and remains untested by a second laboratory. No absolute-level benchmark row is available
either, because the concentration-time courses are figure-only and only the six end-point levels in
section 3's prose table are printed.

## 5. Flags

1. **The second-order unit is per kilogram of food, and the ×10³ convention on that column is only
   inferred.** The header reads `kg × µmol⁻¹ × min⁻¹ × 10³` for steps 5, 9 and 10. For the
   first-order column the ×10³ convention is *confirmed* by the authors' own prose (step 1 printed
   as 6.9 in the table and as 0.0069 min⁻¹ in the Results). **No such confirmation exists for the
   second-order column**, and its printed values (0.0009 to 0.00004) are already tiny, so a reader
   who applies the ×10⁻³ a second time by mistake is out by three decades. Working reading:
   k5 at 160 °C = 0.003 × 10⁻³ = **3×10⁻⁶ kg·µmol⁻¹·min⁻¹ = 3×10⁻³ kg·mmol⁻¹·min⁻¹**. Mark it
   inferred wherever it travels. And note the deeper problem: a per-kilogram-of-food second-order
   constant is not a chemical rate constant, it is a rate constant *times an unknown effective
   volume*, and it cannot be moved to any other matrix.
2. **"Total amino acids" is 72 % protein-bound lysine.** 2112 mg/kg free + 5401 mg/kg bound = 7513
   (mine, exact against the printed total). Every second-order constant (k5, k9, k10) is per unit of
   that pool, and every amino-acid time course in Figure 1 is that pool. Wave B20 exists in this
   repository precisely because bound lysine and free amino acid are not interchangeable reactants;
   this paper lumps them. Any use of k5, k9 or k10 must carry the lump.
3. **No thermocouple and no heat-up correction, in whole nuts.** The authors argue from Demir 2002's
   4-6 min heat-up that the correction is unnecessary, and then assert their own is faster on the
   grounds that they had no thermocouples to install — an argument from the absence of a measurement.
   The earliest sample is at 15 min, so the affected fraction is smaller than in Knol's aqueous
   tubes, but a 4-6 min ramp inside a 15 min point is 30 % of the shortest run and it biases the
   short-time constants downward. **Whole nuts also have an internal gradient**: the surface of a
   hazelnut and its centre are not at the same temperature, and the model is isothermal.
4. **Five constants are marked indeterminate by the authors and seven cells sit on a boundary at
   exactly 0 ± 0.** The named indeterminate ones are **k13, k16, k20, k22, k25** — which includes
   **k16, the largest constant in the paper**, and **k25, the glyoxal sink at 170 °C**, both of which
   this dossier flags as high-value. A zero cell is an estimate that reached its bound, not a
   measurement that the reaction does not occur. **k22 (1-DG → P3) is indeterminate at all three
   temperatures and should not be used at all.** Several HPDs approach or exceed their own estimate:
   k13 @160 (±0.030 on 0.022), k19 @160 (±127 on 134), k18 @160 (±26.1 on 29.4), k15 @170 (±0.58 on
   0.56), k10 @160 (±0.00030 on 0.00027), k20 @160 (±5.6 on 4.7).
5. **The activation-energy table is not on disk and the authors say the barriers are not Arrhenius
   anyway.** Table S4 is the only place Ea and the reparametrised k_b appear, and it is absent from
   the PDF; the body text gives only the range 0-1174 kJ/mol with six zeros. My own three-point
   refits reproduce that pathology (90 to 557 kJ/mol, two of them physically impossible) and are
   recorded in section 3 as `derived_assumption` for that reason only. **No barrier from this paper
   may enter any registry.** This is the single most valuable thing to request (Flags 10).
6. **60 to 85 % of the tracked moles disappear into unmeasured P-pools.** The mass balance is
   printed: 39 %, 29 % and 15 % recovery at 120 min. Every sink constant k20-k26 is fitted to a
   parent's disappearance with nothing on the product side, so read them as loss rates and never as
   named reactions. This is exactly the caveat `parameters_acrylamide.py` carries for `k_acr_dp` and
   `knol2005_extraction.md` states for k6.
7. **Lipid-derived glyoxal and methylglyoxal are in these numbers and the authors say so.** ">56 %
   oil", and Fujioka & Shibamoto's heated olive oil makes 0.61 mg/kg MGO and ~0.5 mg/kg GO. The
   authors bound the contamination at "not ... more than 10 % of their total concentration" by
   comparing against that olive-oil study, and warn that it may inflate **k16 and k6** (formation)
   or **k23 and k25** (degradation). The bound is an argument from a literature comparator, not a
   measurement in this matrix. **The engine's trunk has no lipid lane at all**, so any glyoxal or
   methylglyoxal level compared against this paper inherits an unquantified oil-derived offset.
8. **Glyoxal has a raw-nut background of 1.7 ± 0.6 mg/kg dw, and the model starts from zero.** It is
   the only α-dicarbonyl present before roasting. Any benchmark built on this paper's glyoxal must
   set a non-zero initial condition, and the printed "up to 4 times after 15 min, unchanged
   thereafter" means the roasting-derived increment is only about 5 mg/kg — a small difference
   between two larger numbers, fitted by a sink constant that the authors mark indeterminate at
   170 °C.
9. **Semi-quantitation of three of the seven dicarbonyls.** Only glyoxal, methylglyoxal and
   dimethylglyoxal have their own quinoxaline standards; 3-DG has its own derivatised curve; **1-DG
   and 3,4-DG share ion 235 / 217 and are quantitated against 3-DG's response**, and glucosone (SIM
   251) has no standard listed either — which matters for its reported absence (a null on a
   compound with no calibrant is weaker than a null on a calibrated one, though the m/z is
   monitored). Every constant touching 1-DG (k12, k15, k16, k17, k22) and 3,4-DG (k18, k19) carries
   an unknown multiplicative scale, the same caveat `k_tdg_ddg` and `k_ddg_hmf` already carry.
10. **What to request from the authors**: (i) **Table S4 — the activation energies and k_b at
    160 °C**, without which not one barrier from this laboratory's real-food fit exists on disk;
    (ii) Table S1, the pH and moisture against roasting time, without which every constant here is
    at an unrecorded pH in an unrecorded water activity; (iii) the numeric data behind Figure 1
    (eleven responses × three temperatures × five times); (iv) Table S3, the colour values, which
    would be the first browning response from this laboratory and the only route to testing the
    trunk's browning hold-out against a real food; (v) confirmation of the ×10³ convention on the
    second-order column; (vi) whether glucosone was below a stated limit of detection or simply
    absent from the chromatograms.
11. **Registry gaps against `data/keys/compounds.yml`.** Present: `hmf`; `2_3_butanedione` (this
    paper's DMG). Absent and needed before any benchmark from this paper could be keyed: **glyoxal,
    methylglyoxal, glucosone, 1-deoxyglucosone, 3-deoxyglucosone, 3,4-dideoxyglucosone**, and the
    reactants **sucrose, glucose, fructose** and a lysine key for the bound pool (`reactive_lysine`
    exists as a family id but is not a compound). Six of the eight species whose *levels* this paper
    prints have no registry id. Note also that the repository's own species name for diacetyl is
    `DA` while the registry's is `2_3_butanedione` and this paper's is `DMG`: three names for one
    molecule, and the mapping is not written down anywhere outside this dossier.
12. **What this paper does not contain**: any melanoidin or browning measurement; any pH or a_w
    number on disk; any barrier on disk; any tabulated concentration; any second-order constant in a
    volume basis; any water-activity or moisture series; any variation of the amine (one lumped
    pool); any replicate count for the kinetic runs (the concentration errors are printed as ± but
    n is never stated); any temperature below 150 °C, which is the whole window the trunk was fitted
    in.
