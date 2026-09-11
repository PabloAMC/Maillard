# Zhai et al. 2023 (JAFC) — EXTRACTION (purified TTCA 10 mmol/L in water with and without added L-cysteine at 0.1-10 mol/mol, pH 5.5-8.5, 90-150 C, 20-180 min; browning by A420, volatiles by HS-SPME-GC-MS, ARP and four alpha-dicarbonyls by HPLC, plus five ammonium-sulfide validation models at 120 C — a MECHANISM paper with no rate constant and no activation energy)

### THE MASS-ACTION TEST OF THE STEP THE FIT CANNOT GET RIGHT: adding free cysteine to a TTCA pot delays and suppresses ARP, both deoxypentosones and both short dicarbonyls, which says the engine's one-way `r_ttca_cys` (TTCA -> Cys + pentose) is really an equilibrium — the same conclusion Zhai 2021 reached by pushing on the *xylose* side, now reached by pushing on the *cysteine* side, and the B16 pre-registration's diagnosis was that this exact step runs about ten times too fast.

**Source on disk:** `data/articles/Zhai2023b.pdf` (12 pp., J. Agric. Food Chem. 2023, 71 (39),
14300-14311, doi 10.1021/acs.jafc.3c04166). Read from the `pdftotext -layout` text layer
(`scratchpad/articles/Zhai2023b.txt`, 695 lines), which came through clean for the whole running text,
Materials and Methods, abbreviations, funding and all 44 references. **This paper contains NO table in
the main text.** Every dataset is a figure: Fig. 1 (a, b, c — browning against Cys amount, against
seven temperatures, against seven pH values), Fig. 2 (a-c — A420 time courses at 100/120/140 °C),
Fig. 3 (a-c — the volatile grid at 100/120/140 °C; d — the proposed pathway scheme), Fig. 4 (a-c —
sulfur-, nitrogen- and oxygen-containing totals), Fig. 5 (a-e — ARP, 3-DX, 1-DX, MGO, GO time
courses), Fig. 6 (the five validation models). Under this format's house rule **figures are
`figure_only` and are not typed as numbers**, so what section 3 re-types is the complete set of values
printed in the *running text*, which is all this paper prints outside its figures. The Supporting
Information (Fig. S1, Table S1, Table S2, Fig. S2, and a MGO-GO-Cys flavour table) is **not on disk**.

**⚠ NAMING COLLISION — READ THIS BEFORE USING EITHER FILE.** This paper **already has a dossier**,
written 2026-08-29 in the wave-K6a format under the name **`zhai2023jafc_extraction.md`** (954 lines).
That file is named by *paper identity*; this one is named by *file stem*, as the present brief
requires. **They are the same paper.** Its section 0 also records that the two Zhai 2023 PDFs on disk
are swapped relative to their filenames: `Zhai2023b.pdf` is the **JAFC** paper (this one) and
`Zhai2023.pdf` is the **Food Chemistry** paper (`zhai2023foodchem_extraction.md`). The sibling dossier
digitised Figs. 1b, 3, 5 and 6 and marks every such value `[D]`; this dossier does **not** repeat those
values, because digitised figure readings are `figure_only` here. Use `zhai2023jafc_extraction.md` for
the digitisation work; use this file for the printed record, the methods audit, and the comparison
against Zhai 2021. **Do not treat them as two papers.**

## 0. Identity

| field | value |
|---|---|
| Title | "Regulated Formation of Inhibited Color and Enhanced Flavor Derived from Heated 2-Threityl-Thiazolidine-4-Carboxylic Acid with Additional Cysteine Targeting at Different Degradation Stages" |
| Authors | Yun Zhai, Heping Cui, Khizar Hayat, Teng Li, Xian Wu, **Yuying Fu** (corresponding, webfu@126.com), **Xiaoming Zhang** (corresponding, xmzhang@jiangnan.edu.cn), **Chi-Tang Ho** (corresponding, ctho@sebs.rutgers.edu) |
| Affiliations | Zhejiang Gongshang University (Zhai, Li, Fu); State Key Laboratory of Food Science and Resources, Jiangnan University, Wuxi (Cui, Zhang); Miami University, Ohio (Hayat, Wu); Rutgers University (Ho) |
| Venue | J. Agric. Food Chem. 2023, 71 (39), 14300-14311. Received 20 June 2023, revised 31 August 2023, accepted 13 September 2023, published 25 September 2023 |
| DOI | 10.1021/acs.jafc.3c04166 |
| Paper type | **Mechanism / process-control.** No kinetic model is fitted, **no rate constant and no activation energy appears anywhere in this paper.** Contrast Zhai 2021, which fits three zero-order constants and prints an Ea. |
| Funding | NSFC 32172330; National Key R&D Program of China 2017YFD0400105; Jiangsu Postgraduate Research & Practice Innovation Program KYCX20_1880; National First-Class Discipline Program of Food Science and Technology JUFSTR20180204 |
| Abbreviations as the paper defines them | **TTCA** 2-threityl-thiazolidine-4-carboxylic acid; **ARP** Amadori rearrangement product; **1-DX** 1-deoxy-2,3-pentodiulose; **3-DX** 3-deoxypentose-2-ulose; **GO** glyoxal; **MGO** methylglyoxal; **MRI** Maillard reaction intermediate |
| Direct predecessors by the same group | ref 12 = Zhai et al., JAFC 2019, 67:8632 ("Interaction of added L-cysteine with TTCA ... affecting its Maillard browning") — **the parent of this paper**; ref 13 = Zhai et al., JAFC 2020, 68:10884 (TTCA <-> ARP transformation, `zhai2020_extraction.md`); ref 14 = Zhai et al., JAFC 2021, 69:10648 (TTCA + extra **xylose**, `zhai2021_extraction.md`); ref 37 = Zhai et al., Food Chem. 2023, 404:134420 (`zhai2023foodchem_extraction.md`) |
| Companions on disk | `zhai2023jafc_extraction.md` (**the same paper**, see the collision note above), `zhai2021_extraction.md`, `zhai2020_extraction.md`, `zhai2023foodchem_extraction.md`, `kang2026_extraction.md` and `kang2026_SI_extraction.md` (whose Table S4 the sibling dossier identifies as this family's data re-published) |

## 1. Why it matters

The engine gives TTCA two fates and both are one-way:

| row in `sulfur.py` (lines 550-583) | rate key | order | what this paper bears on |
|---|---|---|---|
| `r_ttca_cys`: TTCA -> Cys + PENT (ring-opening / retro-condensation) | `k_ttca_cys` (fitted on Kang 2026's free-cysteine readings) | 1 | **directly, and adversely** — see below |
| `r_ttca_deg`: TTCA -> 1-deoxypentosone + fragments, no free cysteine released | `k_ttca_deg` (fitted) | 1 | indirectly, through the deoxyosone time courses |

**The B16 diagnosis is the reason to read this paper.** `results/validation/kinetic_core_b16_prereg.md`
records that the shipped lane consumes TTCA about **18x faster** than Zhai 2021 measured at 120-140 °C,
and localises the fault precisely: *"the RING-OPENING step `r_ttca_cys` (TTCA -> Cys + pentose, fitted
on Kang 2026's free-cysteine readings) runs about ten times too fast, and the released pentose is then
eaten by the sugar trunk. The fitted `k_ttca_deg` is not the issue; `k_ttca_cys` is."* Probed on the
shipped B9 lane, TTCA 10 mM at pH 7 after 60 min leaves 1.58 / 0.05 / 0.00 mM at 100 / 120 / 140 °C
against Zhai 2021's measured 8.71 / 6.07 / 4.46 mM.

**What this paper adds to that diagnosis is the missing half of a mass-action argument.** Zhai 2021
pushed on the *xylose* side and found that extra xylose **accelerates** TTCA degradation, by trapping
the released cysteine and pulling the ring-opening forward. This paper pushes on the *cysteine* side
and finds the mirror image: added free cysteine **delays and suppresses** every downstream marker —
"At 120 or higher temperature (140 °C), 3-DX was undetectable before the reaction period of 40-60 min
in the TTCA-Cys model, confirming that the degradation pathway of TTCA to form downstream products was
blocked", and "ARP was not observed at the early stage, which could confirm the phenomenon of Cys
blocking the reaction process through its action on the targeted TTCA **to make it convert back to Cys
and Xyl**." Two perturbations, opposite directions, same conclusion: **`r_ttca_cys` is an equilibrium,
and the engine has no reverse for it.** A one-way ring-opening fitted against a free-cysteine reading
in a pot where the reverse is active will come out too fast, which is exactly the B16 residual. This
is a **structural** finding — the paper supplies no equilibrium constant and no reverse rate — but it
is the first evidence in the corpus that the fault B16 localised has a named mechanism rather than a
bad fit.

**Second, the two furan-to-thiol edges are validated qualitatively in the same pot type.** Figure 6's
validation models charge ammonium sulfide (10 mmol/L) against furfural, furan, 2-methylfuran and
**4-hydroxy-5-methyl-3(2H)-furanone** — which is **norfuraneol**, the engine's `NF` — at pH 7.0 and
**120 °C**, and report that the four furanoids "were rapidly consumed in their respective thermal
reaction models and were all exhausted after 20 min", giving respectively 2-furfurylthiol,
2-methylthiophene, thiophene and **2-methyl-3-furanthiol**. Two of those four are engine reactions:

- `r_fur_fft`: FUR + H2S -> FFT (`k_fur_fft`, `neutral_h2s`) — panel (a);
- `r_nf_mft`: NF + H2S -> MFT (`k_nf_mft`, `neutral_h2s`) — panel (b).

The other two products, thiophene and 2-methylthiophene, are compounds the module **declares out of
scope**. The paper's mechanistic claim for all four is the same: "the oxygen atoms on the ring or the
side chain of the ring in the furans were replaced by sulfur atoms."

**Third, it says how small the volatile lane is.** The total quantified volatile pool at 180 min is
**87.262 µg/L** at 140 °C with added cysteine, against a 10 mmol/L TTCA charge. On a rough molar basis
that is of order **0.008 mol %** of the charged intermediate (section 3, arithmetic 4, mine). The
aroma channel is a trace branch off a pot whose carbon overwhelmingly goes to browning — a useful
calibration for how much flux the sulfur lane's product steps should be expected to carry.

What this paper does NOT give: any rate constant; any activation energy; any TTCA concentration at any
time (it measures TTCA's *products*, never TTCA itself); any thiol dosing experiment; any disulfide
measurement; any measurement of any thiol **consumption** channel; any pH-resolved volatile data (the
pH ladder is browning-only); any error bar on the printed volatile totals; and any of its Supporting
Information.

## 2. Methods as they matter to a model

- **TTCA synthesis and purification.** "The solution of Cys and Xyl (**0.0827 mol/L**) with the same
  molar ratio was prepared using deionized water and adjusted to **pH 7.4 ± 0.01** by NaOH solution
  (2 and 6 mol/L). The solution was then heated under **90 °C for 40 min** and subsequently cooled by
  an ice bath." Purified through a **Dowex 50WX4** ion-exchange resin in H+ form (200-400 mesh) and by
  semi-preparative RP-HPLC on an **Xbridge amide** column (4.6 mm × 150 mm, 3.5 µm, Waters); procedure
  as in ref 13. Identity confirmed by **UPLC-ESI-MS** (Waters Synapt MALDI Q-TOF) and NMR (Bruker DRX
  400) — spectra in Fig. S1 and Table S1, **not on disk**. **No purity figure is printed in this
  paper** (Zhai 2021 states 98 %).
- **The reaction pots.** "TTCA reaction models (**10 mmol/L**) with or without the addition of
  different amounts of Cys (**0.1, 0.5, 1, 2, and 10 mol/mol TTCA**) were prepared under different pH
  values (**5.5, 6, 6.5, 7, 7.5, 8, and 8.5**) and heated at different temperatures (**90, 100, 110,
  120, 130, 140, and 150 °C**) for varying times (**20, 40, 60, 80, 100, 120, 140, 160, and
  180 min**). The reaction was performed in appropriate **pressure-resistant bottles** equipped with a
  collector-type **magnetic stirrer (DF-101S)**. The reaction was terminated in an **ice bath**."
  **Aqueous, unbuffered** — the pH is set by adjustment, as in Zhai 2021 (NaOH/HCl); no buffer is
  named anywhere in the reaction section.
- **The three grids are not the same grid.** The full seven-temperature and seven-pH ladders exist
  **for browning only** (Fig. 1b at 120 min, pH 7, Cys 1:1; Fig. 1c at 120 °C, 120 min, Cys 1:1). All
  volatile, ARP and dicarbonyl data are at **100, 120 and 140 °C only** and at **pH 7.0, Cys : TTCA
  1 : 1**. This asymmetry is the single most important thing to know before quoting "seven
  temperatures" about anything in this paper (Flags 2).
- **Browning.** A420 on a Shimadzu 2100 UV-vis; samples "diluted to an appropriate concentration to
  make the measured value fall into the appropriate range" — **the dilution factors are not printed**,
  so the A420 values are on an unstated basis.
- **ARP.** HPLC-ELSD, **external standards**, XBridge BEH amide (3.5 µm, 4.6 × 150 mm), linear
  gradient of 10 mM ammonium formate (A, pH 6) and 100 % acetonitrile (B), 1 mL/min, 10 µL injection,
  column 25 °C; ELSD drift tube 45 °C, N2 1.5 L/min. **"The linear gradient elution process is not
  shown here."**
- **alpha-Dicarbonyls (3-DX, 1-DX, MGO, GO).** Derivatised with **OPD** (2 g/100 mL) plus **DTPA**
  (11 mmol/L) in HEPES at pH 7; sample and reagent mixed 1 : 1 and incubated **in the dark at 25 °C
  for 4 h**; the quinoxaline derivatives quantified against **external standards**; HPLC-DAD, Waters
  SunFire C18 (5 µm, 4.6 × 150 mm). **"Linear gradient elution process is not presented here."**
- **Volatiles — HS-SPME-GC-MS.** **3 g** of sample plus **3 µL** of internal standard
  (**1,2-dichlorobenzene, 0.018 µg/µL in methanol**) in a 20 mL headspace vial with a PTFE/BYTL
  septum; **CAR/PDMS/DVB 75 µm** fibre exposed to a depth of **2.5 cm**; extraction **60 °C for
  20 min** in a thermostatic water bath; desorption **10 min at 250 °C**; Agilent **7890B** GC with
  **5977B** MSD; **DB-Wax** 30 m × 0.25 mm × 0.25 µm. **The oven programme is not printed** — "set and
  applied according to our previous research", i.e. deferred to Zhai 2021. Identification by **NIST 17**,
  Kovats retention indices, **WILEY 07** and literature; quantification by **calibration curves
  (Table S2, not on disk)** with x = concentration (µg/L) and y = the peak-area ratio of standard to
  internal standard.
- **The validation models.** "Reaction models of ammonium sulfide-furfural, ammonium sulfide-furan,
  ammonium sulfide-2-methylfuran, ammonium sulfide-4-hydroxy-5-methyl-3(2H)furanone, and ammonium
  sulfide-MGO-GO were prepared. Compounds ... were added to the aqueous **ammonium sulfide solution
  (10 mmol/L)** with the **same molar ratio** ... adjusted to **pH 7.0 ± 0.1** and subsequently treated
  at **120 °C** for a defined reaction time (**10, 20, 30, and 40 min**)", in pressure-resistant
  bottles, ice-quenched, analysed by GC-MS.
- **Statistics.** "All measurements were conducted in **triplicate**, and results were presented as
  mean ± standard deviation (SD). Significant differences (**p < 0.05**) were calculated by ANOVA
  using SPSS 21 ... All the experiments in the research process were **performed three times**."
  **No SD is attached to any number printed in the running text.**
- **Temperature window.** 90-150 °C covers the sulfur lane's own 95-145 °C window exactly — the widest
  such coverage in the TTCA corpus — but only for A420.

## 3. Tables re-typed

**There is no table in the main text of this paper.** Tables S1 and S2 are in the Supporting
Information, which is not on disk. What follows is therefore the complete set of quantities printed in
the **running text**; every figure-borne value is left in its figure per the house rule.

### 3.1 Quantitative results printed in the running text

| quantity | value | conditions | where |
|---|---|---|---|
| **total volatile flavour, TTCA-Cys model, at the "later reaction period"** | **29.565 / 63.249 / 87.262** µg/L | 100 / 120 / 140 °C, pH 7.0, Cys : TTCA 1 : 1 | Results, "Additional Cys as the Inhibitor..." |
| **total volatile flavour, TTCA model (no added Cys)** | **13.749 / 39.168 / 63.035** µg/L | same temperatures, pH 7.0 | same paragraph |
| number of **types** of sulfur-containing volatiles at the later period, **100 °C** | **27** (TTCA-Cys) against **19** (TTCA) | 100 °C, pH 7.0 | same paragraph |
| number of pyrazine species detected in the TTCA-Cys model | **four kinds** | (temperature not restated — see Flags 4) | same paragraph |
| **total pyrazine content at 180 min** | **1.203** µg/L (TTCA-Cys) against **0.739** µg/L (TTCA) | 180 min, pH 7.0 | same paragraph |
| optimal added-Cys loading for colour inhibition | **1 mol/mol TTCA** ("the Cys addition of **2 mol/mol** TTCA presented the best color-inhibiting effect, but the browning index ... did not show a significant downward trend compared to that with the Cys addition of **1 mol/mol** TTCA. Therefore, **1 mol/mol TTCA was determined as the optimal amount**") | 120 °C, pH 7, 120 min | Results, first section |
| threshold above which Cys blocks melanoidin formation | **≥ 0.5 mol/mol TTCA** | 120 °C, pH 7, 120 min | same |
| the loading at which colour inhibition **reverses** | **10 mol/mol TTCA** ("the color inhibition effect showed a downward trend instead") | same | same |
| the loading at which added Cys **increases** colour slightly | **1 : 10 Cys : TTCA** ("revealed a little bit increase compared to the TTCA model") | same | same |
| pH range over which TTCA browning barely moves | **5.5-7.5** ("a small fluctuation trend") | 120 °C, 120 min | Results, first section |
| pH range over which TTCA browning rises significantly | **7.5-8.5** | same | same |
| the four validation furanoids at 120 °C | "rapidly consumed ... and were **all exhausted after 20 min**" | (NH4)2S 10 mmol/L, equimolar substrate, pH 7.0, 120 °C | Results, "Model Reaction Systems..." |
| 3-DX detectability in the TTCA-Cys model | **undetectable before 40-60 min** | 120 and 140 °C, pH 7.0, 1 : 1 | Results, "Dynamic Formation..." |
| ARP in the TTCA-Cys model at 100 °C | "slightly higher than that in the TTCA reaction model during the **initial 60 min**" | 100 °C | same |
| TTCA synthesis charge | Cys + Xyl **0.0827 mol/L**, equimolar, pH **7.4 ± 0.01**, **90 °C / 40 min** | — | Methods |
| reaction charge | TTCA **10 mmol/L** | — | Methods |
| validation charge | (NH4)2S **10 mmol/L**, equimolar substrate | pH 7.0 ± 0.1 | Methods |
| internal standard | 1,2-dichlorobenzene, **0.018 µg/µL**, **3 µL** into **3 g** of sample | — | Methods |

### 3.2 Products named as increasing with added cysteine (presence and direction only, no numbers printed)

At 100 °C, "the contents of **2-furanthiol, 2-methyl-3-furanthiol, 2-acetylthiazole, and thiophene**
exhibited a substantial increase"; **thiophene[3,2-b]thiophene** and
**2-methylthiophene[2,3-b]thiophene** "were not found in the TTCA reaction model, but additional Cys
could significantly promote their generation." Furfural is named as the dominant oxygen heterocycle
("furfural exhibited the highest concentration") and its **suppression** by added Cys is the paper's
central flavour observation: "As the temperature increased to 120 or 140 °C, only a small amount of
furans could be detected at 40 or 80 min and not even detected as the time extended to 180 min ...
these furans were captured by a large number of active H2S generated from Cys degradation."

### 3.3 The validation-model product pairs, as the text states them

| substrate | product named | engine reaction, if any |
|---|---|---|
| furfural | **2-furfurylthiol** | **`r_fur_fft`** (FUR + H2S -> FFT) |
| **4-hydroxy-5-methyl-3(2H)-furanone** (= norfuraneol, engine `NF`) | **2-methyl-3-furanthiol** | **`r_nf_mft`** (NF + H2S -> MFT) |
| 2-methylfuran | 2-methylthiophene | none — **out of the module's declared scope** |
| furan | thiophene | none — **out of scope** |
| MGO + GO | "a large amount of **thiazoles** and **thiophenes** and a small amount of nitrogen-containing flavor substances" | partially `r_cys_actz` (Cys + MGO -> 2-acetylthiazole) |

### Arithmetic on the printed numbers (all mine)

**1. What added cysteine does to the total volatile pool, and how it fades with temperature.**
TTCA-Cys / TTCA at 180 min: 29.565/13.749 = **2.15** (100 °C); 63.249/39.168 = **1.61** (120 °C);
87.262/63.035 = **1.38** (140 °C). **The cysteine boost shrinks monotonically as temperature rises** —
which is the arithmetic form of the paper's own claim that at high temperature the TTCA pot generates
its own sulfide fast enough that an external cysteine supply matters less.

**2. Temperature response of each arm.** TTCA alone: ×2.85 over 100 -> 120 °C, ×1.61 over
120 -> 140 °C, **×4.58 overall**. TTCA + Cys: ×2.14, ×1.38, **×2.95 overall**. **The added-cysteine
arm is the flatter of the two in temperature**, so the two curves converge; extrapolating either fold
outside 100-140 °C is not supported by anything here.

**3. Pyrazines.** 1.203/0.739 = **×1.63** with added cysteine. As a share of the total volatile pool
at 180 min (using the 100 °C totals, on the reading that the pyrazine sentence sits in the 100 °C
paragraph — Flags 4): **4.1 %** (TTCA-Cys) against **5.4 %** (TTCA). So added cysteine raises the
pyrazines in absolute terms while *lowering* their share, because it raises the sulfur compounds more.

**4. The whole volatile pool is a trace branch (order of magnitude, mine).** TTCA is
cysteine + xylose − water = 121.16 + 150.13 − 18.02 = **253.27 g/mol**, so 10 mmol/L is **2.53 g/L**.
Taking a representative volatile molar mass of ~114 g/mol (both MFT and FFT are C5H6OS = 114.17), the
largest printed total, 87.262 µg/L, is ~**0.77 µmol/L** against 10 000 µmol/L of charged TTCA, i.e.
**~0.008 mol %**. Even at ×10 for the molar-mass assumption and the compounds outside the calibrated
set, the quantified aroma pool is **well under a tenth of a mole percent** of the intermediate. The
carbon is in the browning. **Assumptions declared: one representative molar mass; the sum covers only
the calibrated compounds; HS-SPME quantification against one internal standard.**

**5. A bound on the furanoid consumption rate at 120 °C (mine, and weak).** "All exhausted after
20 min" on a 10/20/30/40 min sampling grid means only *below the detection limit at the 20 min
sample*. If "exhausted" is taken as ≥ 95 % consumed, the pseudo-first-order constant is
≥ ln(20)/20 min = **0.15 min⁻¹** against 10 mmol/L sulfide; at ≥ 99 % it is ≥ 0.23 min⁻¹. **This is
not a measurement**: the detection limit is not stated, the disappearance includes every non-thiol
route, the sulfide activity of (NH4)2S at pH 7 and 120 °C is not stated, and the substrate is at
10 mmol/L rather than at the µg/L levels of the Maillard pot. Recorded as a floor, class
`derived_assumption`. For scale, Yaghmur 2005 measures furfural in water with cysteine at 65 °C with a
half-life of 5.5 h; this is the same direction, 55 °C hotter, with a sulfide charge.

## 4. Kinetic numbers the repository can use

**There is no rate constant and no activation energy in this paper.** The usable content is (a) three
pairs of printed volatile totals and their ratios, (b) a set of directional and threshold statements
about the added-cysteine perturbation, and (c) the qualitative validation of two engine edges.

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Keyed: `2_methyl_3_furanthiol`,
`2_furfurylthiol`, `furfural`, `hydrogen_sulfide`. **Not keyed:** TTCA, the xylose-cysteine ARP,
cysteine, xylose, 3-DX, 1-DX, methylglyoxal, glyoxal, 2-acetylthiazole, 2-furanthiol, thiophene,
2-methylthiophene, thieno-thiophenes, 4-hydroxy-5-methyl-3(2H)-furanone (norfuraneol), or any
pyrazine.

Every row below shares: **purified TTCA 10 mmol/L in unbuffered water**, ± L-cysteine, ice-quenched in
pressure-resistant stirred bottles, triplicate, ANOVA at p < 0.05. Unless stated otherwise the
condition is **pH 7.0, Cys : TTCA = 1 : 1**.

| step / observable | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| total volatiles, TTCA-Cys | sum over calibrated compounds at the later period | **29.565 / 63.249 / 87.262** | µg/L | 100 / 120 / 140 °C, pH 7, 1 : 1 | **nothing is fitted in this paper** | Results text, p. 14306 | **level_only** (an HS-SPME sum against one internal standard; see Flags 3 on its five-figure precision) |
| total volatiles, TTCA alone | same | **13.749 / 39.168 / 63.035** | µg/L | same temperatures | — | same | **level_only** |
| effect of added cysteine on the volatile pool | ratio TTCA-Cys / TTCA | **2.15 / 1.61 / 1.38** at 100 / 120 / 140 °C | — | as above | — | derived from the printed totals (mine) | **within_study_ratio** — the strongest quantitative statement in the paper |
| temperature response, TTCA alone | fold 100 -> 120 -> 140 °C | 2.85 / 1.61 (4.58 overall) | — | as above | — | derived (mine) | **within_study_ratio** |
| temperature response, TTCA + Cys | fold 100 -> 120 -> 140 °C | 2.14 / 1.38 (2.95 overall) | — | as above | — | derived (mine) | **within_study_ratio** |
| total pyrazines at 180 min | TTCA-Cys vs TTCA | **1.203** vs **0.739** (ratio 1.63) | µg/L | pH 7, 1 : 1, 180 min, temperature not restated (Flags 4) | — | Results text | **level_only** / the ratio **within_study_ratio** |
| sulfur-compound diversity | number of species detected | **27** vs **19** | count | 100 °C, later period | — | Results text | **level_only** — a count of detections, **not an amount** |
| optimal colour-inhibiting Cys loading | Cys : TTCA | **1 : 1** (2 : 1 best but not significantly better; ≥ 0.5 : 1 blocks melanoidins; 10 : 1 reverses; 0.1 : 1 slightly increases colour) | mol/mol | 120 °C, pH 7, 120 min | — | Results, Fig. 1a discussion | **threshold** (a dose-response threshold set stated in words) |
| pH response of TTCA browning | flat 5.5-7.5, rises 7.5-8.5 | — | — | 120 °C, 120 min | — | Results, Fig. 1c discussion | **level_only** (direction only; **all values are in the figure**) |
| **`r_ttca_cys` is reversible** | added free cysteine delays ARP, 3-DX, 1-DX, MGO and GO and "make[s] it convert back to Cys and Xyl" | — | — | 100-140 °C, pH 7, 1 : 1 | — | Results, "Dynamic Formation...", Fig. 5 | **threshold** (a presence/absence and delay result) — **structural, and the most useful line in the paper for this engine** |
| 3-DX suppression | undetectable before **40-60 min** in the TTCA-Cys model | min | | 120 and 140 °C | — | Results text | **threshold** |
| **`r_fur_fft` validated** | furfural + (NH4)2S gives 2-furfurylthiol; furfural exhausted by 20 min | — | — | 10 mmol/L each, pH 7.0, **120 °C** | — | Results, Fig. 6a | **threshold** |
| **`r_nf_mft` validated** | 4-hydroxy-5-methyl-3(2H)-furanone (norfuraneol) + (NH4)2S gives 2-methyl-3-furanthiol; substrate exhausted by 20 min | — | — | as above | — | Results, Fig. 6b | **threshold** |
| furanoid consumption floor | pseudo-first-order k | **≥ 0.15** (at ≥95 % consumed) | min⁻¹ | (NH4)2S 10 mmol/L, pH 7.0, 120 °C | first order in the furanoid, assumed | derived from "exhausted after 20 min" (mine) | **derived_assumption** — a floor with an undefined detection limit; **do not fit against it** |
| volatile pool as a share of charged TTCA | ~0.008 | mol % | 140 °C, 180 min, pH 7, 1 : 1 | — | derived (mine, one molar-mass assumption) | **derived_assumption** |
| Figs. 1a-c, 2a-c, 3a-d, 4a-c, 5a-e, 6 — every browning value, every volatile concentration, every ARP / 3-DX / 1-DX / MGO / GO time course, every validation-model bar | — | — | — | — | Figs. 1-6 | **figure_only** — digitised in `zhai2023jafc_extraction.md` and marked `[D]` there; **not typed as numbers here** |

### What is new here, against `zhai2021_extraction.md`

Read side by side, the two papers are the same pot with the opposite perturbation, and only one of
them has kinetics.

| axis | **Zhai 2021** (`zhai2021_extraction.md`) | **this paper (2023 JAFC)** |
|---|---|---|
| perturbation | extra **xylose** (10 mmol/L) | extra **cysteine** (0.1 / 0.5 / 1 / 2 / 10 mol/mol) |
| direction of the effect | **accelerates** TTCA degradation and **deepens** browning, by trapping released cysteine | **retards** the downstream path and **lightens** browning, by pushing the ring-closure back |
| TTCA measured directly? | **yes** — TTCA remaining against time | **no.** TTCA itself is never quantified in this paper; only its products |
| rate constants | **three**, zero-order, c = c0 − kt at pH 7: 0.0271 / 0.0651 / 0.0813 mmol L⁻¹ min⁻¹ at 100 / 120 / 140 °C, R² 0.9516 / 0.9949 / 0.9802 | **none** |
| activation energy | printed **80.99 kJ/mol**, which the 2021 dossier flags as **not reproducible** from the three printed k (re-derivation gives ~35 kJ/mol) | **none** |
| temperature grid | 100 / 120 / 140 °C | **90-150 °C in 10 °C steps — but for browning only**; everything else is 100 / 120 / 140 °C |
| pH grid | 5.5 / 6 / 7 / 8 | **5.5 / 6 / 6.5 / 7 / 7.5 / 8 / 8.5 — browning only, at 120 °C** |
| time grid | 20-120 min | **20-180 min** |
| volatiles | not quantified | **quantified**, 39 compounds by HS-SPME-GC-MS against 1,2-dichlorobenzene, at three temperatures |
| dicarbonyls | 3-DX, 1-DX, MGO, GO against time at **100 °C** and four pH values | the same four **plus ARP**, against time, at **three temperatures**, ± cysteine |
| isotope work | **yes** — ¹³C₅-xylose tracer, and an LC-MS/MS fragment table | **none** |
| mechanism validation | none | **yes** — five ammonium-sulfide models at 120 °C, including furfural -> FFT and norfuraneol -> MFT |
| what it gives the fit | **three measured rates on `k_ttca_deg`'s species** — declared FIT rows in B16/B17/B25 (TTCA remaining after 60 min: 8.71 / 6.07 / 4.46 mmol/L) | **no row.** Its contribution is **structural**: the reverse of `r_ttca_cys` exists and matters |

**In one sentence:** Zhai 2021 gives the fit its numbers on this step; this paper gives the fit its
diagnosis of why those numbers are missed.

### Can any of this be put on the same basis as the sulfur lane's constants?

**(a) The reversibility finding — yes, as a structural change, not a number.** The engine's
`r_ttca_cys` has no reverse. Adding one would need an equilibrium constant, and **this paper supplies
none** — it supplies only the sign of the effect and the delay times. The honest use is to record the
reverse as a declared missing edge against the B16 residual, in the way `THIOL_CHANNELS` records the
thioether channel's reversibility with Stack 2018's K, and to say plainly that the K for this one does
not exist in the corpus.

**(b) The two validated edges — yes, as confirmation, not as rates.** `r_fur_fft` and `r_nf_mft` are
demonstrated in isolation at 120 °C, which is inside the module's window and is a cleaner
demonstration than any in-pot inference. But the loading is 10 mmol/L against a Maillard pot's µg/L,
the sulfide source is ammonium sulfide rather than cysteine-derived H2S, and the branch fractions are
figure-borne. Confirmation of topology; no constant.

**(c) The volatile totals — usable only as within-study ratios.** The 2.15 / 1.61 / 1.38 boost and the
2.85 / 1.61 temperature folds are dimensionless and internally consistent. The absolute µg/L values
depend on calibration curves that are in an SI not on disk, on a single internal standard, and on an
HS-SPME step whose partition behaviour differs per compound. **Never score an absolute level from this
paper.**

**(d) What cannot be used at all.** The browning ladders (all figure-borne, and on an unstated dilution
basis), the pH ladders (figure-borne, browning only), the ARP and dicarbonyl time courses
(figure-borne), and every validation-model amount (figure-borne). Nothing here bears on any entry in
`NO_MEASURED_RATE` — `k_oligomer` remains the sole entry and this paper adds nothing to it — and
nothing here bears on **any** of the six channels in `THIOL_CHANNELS`: there is no thiol dosing, no
disulfide measurement, no thiol mass balance and no thiol-consumption experiment anywhere in the paper.

## 5. Flags

1. **This paper already has a dossier under a different name.** `zhai2023jafc_extraction.md`
   (954 lines, wave K6a, 2026-08-29) is the same paper, and it additionally records that the two Zhai
   2023 PDFs are **swapped relative to their file stems** (`Zhai2023b.pdf` is the JAFC paper;
   `Zhai2023.pdf` is the Food Chemistry paper). Any pipeline keyed on the filename has them the wrong
   way round. **Reconcile the two dossiers before either is cited as an independent source**, and do
   not let the same paper enter an objective twice under two names.
2. **"Seven temperatures" applies to browning only.** The 90-150 °C ladder (Fig. 1b) and the pH
   5.5-8.5 ladder (Fig. 1c) are **A420 measurements**. Every volatile, ARP and dicarbonyl measurement
   is at **100, 120 and 140 °C, pH 7.0, Cys : TTCA 1 : 1**. There is **no seven-rung MFT or FFT ladder
   in this paper**.
3. **The volatile totals are printed to five significant figures with no error bar.** 29.565, 63.249,
   87.262, 13.749, 39.168, 63.035, 1.203 and 0.739 µg/L are sums over a calibrated compound set,
   obtained by HS-SPME against a single internal standard, in triplicate — and **not one of them
   carries an SD**, although the Methods promise "mean ± standard deviation". Three decimals on a
   headspace sum is precision the method cannot support. Use the ratios; treat the absolutes as
   one- or two-figure quantities.
4. **The pyrazine sentence's temperature is ambiguous.** "It is worth noting that four kinds of
   pyrazines could be detected in the TTCA-Cys reaction model, the total content of which increased to
   1.203 µg/L compared with that of 0.739 µg/L in the TTCA system at the final reaction stage
   (180 min)" sits in the paragraph that opens "At 100 °C, the types of sulfur-containing flavor
   compounds increased to 27...". **My reading is that it is the 100 °C figure**, and my share
   calculation in section 3 says so explicitly, but the sentence does not restate the temperature.
5. **Two chromatographic methods are unreproducible as printed.** "The linear gradient elution process
   is not shown here" appears for the ARP HPLC-ELSD method, and "Linear gradient elution process is not
   presented here" for the dicarbonyl HPLC-DAD method. The **GC oven programme** is likewise not
   printed and is deferred to ref 14 (Zhai 2021). Three of the paper's four quantitative methods cannot
   be reproduced from this paper alone.
6. **The retention-index alkane set is stated inconsistently.** Materials says "The n-alkanes
   (**C6-C27**) for retention indices"; Methods says "The calculation of the retention index (RIX) of
   target compounds was performed through **C7-C30** n-alkane standards". One of the two is wrong. No
   retention index is printed in the paper, so nothing downstream depends on it, but it is a sign of
   the level of proof-reading.
7. **TTCA itself is never measured.** The paper infers what happens to TTCA entirely from ARP, 3-DX,
   1-DX, MGO, GO and the volatiles. Anyone hoping this paper would supply a second temperature series
   on `k_ttca_deg` or `k_ttca_cys` — the very thing B16 needs — will not find one. **Zhai 2021 remains
   the only direct TTCA time series in the corpus.**
8. **No purity figure and no equilibrium constant.** The TTCA purity is deferred to ref 13 (Zhai 2021
   states 98 %) and the reverse reaction that this paper's whole argument rests on — cysteine plus
   pentose closing back to TTCA — is **never quantified**: no K, no reverse rate, no measurement of
   TTCA re-formation. The finding is a sign, not a size.
9. **The validation models are three orders of magnitude above the pot.** 10 mmol/L furanoid against
   10 mmol/L ammonium sulfide, versus a Maillard pot whose entire volatile pool is under 0.1 mg/L.
   Ammonium sulfide is also a sulfide *source* whose free-H2S activity at pH 7 and 120 °C the paper
   does not state, so even the *order* in sulfide is unconstrained. **Read Fig. 6 as topology.**
10. **Two of the four validation products are out of the module's scope.** Thiophene and
    2-methylthiophene, from furan and 2-methylfuran, have no species in `species_sulfur.py`, and the
    module declares thiophenes out of scope (see the B8 panel's `zhai_13C_exogenous_carbon_threshold`
    note). This paper says the O -> S substitution runs on **four** furanoids of which the engine
    models two — so the engine's furan-to-thiol chemistry is a subset of a broader class, and the
    sulfide that goes to the thiophenes is, in the engine, unaccounted.
11. **Nothing here touches any thiol-consumption channel.** No thiol is dosed, no disulfide is
    measured, no thiol mass balance is attempted, and the paper reports thiols only as they rise. The
    lane's central defect — the thiol sink, INTRODUCTION section 7 — gets **no evidence** from this
    paper in either direction. That is worth recording as a negative so nobody looks again.
12. **Unbuffered pH, set once.** The reaction solutions are aqueous and pH-adjusted; no buffer is named
    in the reaction section, and pH is not re-measured after heating. In a pot generating H2S, NH3, CO2
    and carboxylic acids from 10 mmol/L of an amino-acid-derived intermediate, the drift is not
    negligible, and the pH ladder's labels are therefore initial values.
13. **What this paper does not contain**: any rate constant; any activation energy; any TTCA
    concentration; any equilibrium constant for the ring-opening; any SD on any printed number; any
    volatile data outside 100/120/140 °C; any volatile data at any pH other than 7.0; any thiol
    consumption measurement; any isotope tracing; any of its Supporting Information (Fig. S1,
    Table S1, Table S2, Fig. S2, the MGO-GO-Cys flavour table).
14. **What to request from the authors**: (i) **TTCA remaining against time in the added-cysteine
    arm** — the one measurement that would turn this paper's structural finding into a number for
    `k_ttca_cys` and its missing reverse; (ii) the equilibrium position of TTCA <-> Cys + xylose at
    100-140 °C, or any measurement of TTCA re-formation from a Cys + Xyl charge at reaction
    temperature; (iii) the Supporting Information, in particular **Table S2's calibration curves**,
    without which the µg/L totals cannot be checked; (iv) the numeric data behind Figs. 3 and 5 (the
    39-compound × 4-time × 3-temperature × 2-system grid, and the five ARP/dicarbonyl series); (v) the
    two omitted HPLC gradients and the GC oven programme; (vi) the SDs the Methods promise; (vii) the
    detection limits behind "exhausted after 20 min" in Fig. 6.
15. **Registry gaps against `data/keys/compounds.yml`**: `2_methyl_3_furanthiol`, `2_furfurylthiol`,
    `furfural` and `hydrogen_sulfide` are present. **Absent: TTCA itself** — the subject of the paper
    and a species the engine carries (`species_sulfur.py`) — together with the xylose-cysteine **ARP**,
    **cysteine**, **xylose**, **3-DX**, **1-DX**, **methylglyoxal**, **glyoxal**, **2-acetylthiazole**
    (an engine species, `ACTZ`), **2-furanthiol**, **norfuraneol** (engine `NF`), **thiophene**,
    **2-methylthiophene**, the two thieno-thiophenes, and every pyrazine. A benchmark row built on
    this paper would need at least TTCA, ARP and norfuraneol keyed.
