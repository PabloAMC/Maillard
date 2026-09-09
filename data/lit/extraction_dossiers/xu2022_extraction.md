# Xu 2022 — EXTRACTION (pea protein isolate 5 g/L in 0.01 M potassium phosphate pH 7.6, four typical cooked-beef aroma compounds SEALED IN THE VIAL AND HEATED WITH THE PROTEIN at 80/90/100/120 C for 10 min or at 100 C for 5-20 min, cooled on ice, equilibrated 12 h at 25 C, read by SPME-GC-MS; Klotz constants, Stern-Volmer quenching and van 't Hoff thermodynamics for 2-methylpyrazine before and after a 100 C / 10 min treatment)

### THE 120 °C PEA EXPERIMENT — AND IT IS DOSE-ADDED-PRE-COOK. Wave B26's record says of the Bi 2022 pea rows: "37 C is an in-mouth temperature, not a process one: nothing here licenses a pea binding constant at 90 or 140 C." This paper is by **the same laboratory, on the same protein, in the same buffer, with Bi Shuang as second author and Wu Jihong as corresponding author** — and it takes pea protein isolate to **80, 90, 100 and 120 °C for 10 min with the aroma compounds sealed in the vial**, then cools and reads the binding at 25 °C. It prints a **pea Klotz constant before and after a 100 °C / 10 min treatment: nK rises from 13 065.25 ± 546.63 to 22 876.70 ± 985.48 L/mol, a 1.75x increase (Table 3, p. 34)** — heating pea protein *increases* its binding of 2-methylpyrazine. And because it prints **both a bound fraction and a protein loading**, its binding rates convert to the registry's own `K_g` form **without needing a molar mass at all**, the obstruction that made Bi 2022's Klotz K unusable in L/g: **2-methylpyrazine on pea = 4.27e-2 L/g native, 6.79e-2 after 100 °C, ~1.97e-1 after 120 °C (mine)**, against the shipped `kg_hexanal_pea` of 2.537e-1. **But the measurement temperature is still 25 °C** — the 120 °C is the cook, not the equilibrium — and **the paper itself says the two highest-binding compounds are bound COVALENTLY** ("由于共价结合是不可逆的" — because covalent binding is irreversible), which quarantines them from `REVERSIBLE_BINDING` on the authors' own testimony.

**Source on disk:** `data/articles/Xu2022.pdf` (8 pp., 食品科学 / *Food Science* (China) 2022, 43(19): 28-35).
**Language: the body is in CHINESE; the title, author list, affiliations, abstract, keywords, table captions,
figure captions and reference list are in ENGLISH.** All three table captions are bilingual. Read from the
`pdftotext -layout` text layer, which carried the Chinese cleanly except for the embedded equation glyphs;
**Tables 1 (p. 32), 2 (p. 34) and 3 (p. 34) were all verified cell-by-cell against the rendered page images**
and are re-typed in full below with the Chinese terms translated. **What is NOT legible or not present**:
Equation (1) rendered as mojibake in the text layer and was recovered from the page image (it is the standard
headspace-depletion form and is given below); the two-character sensory descriptors in Table 1 are Chinese only
and are translated here with the original characters retained. **Figures 1-5 are images and carry no printed
numbers**: Fig. 1 (binding rate against heat-treatment temperature, four compounds × five conditions), Fig. 2
(binding rate against heat-treatment time, four compounds × five conditions), Fig. 3 (surface hydrophobicity,
panels A and B), Fig. 4 (particle size and zeta potential, panels A-D) and Fig. 5 (the Klotz double-reciprocal
plot, two lines, with R^2 = 0.979 and 0.978 printed inside the panel). **The binding rates are therefore almost
entirely figure-only: only FIVE numeric binding rates appear anywhere in the text, plus the ">98 %" bound.**
There is **no supplementary material**. Repo status before this dossier: Xu 2022 has **no extraction dossier**
and is **not cited** in `src/kinetic_core/parameters_matrix.py`, `src/kinetic_core/matrix_sites.py` or
`data/species/protein_matrices.yml`. **Note the name collision on disk**: `xu2010_extraction.md` already
exists and is a different paper.

## 0. Identity

| field | value |
|---|---|
| Title (English, as printed) | "Effect of Heat Treatment on the Binding of Pea Protein Isolate to Typical Beef Aroma Substances and Their Interaction Mechanism" |
| Title (Chinese) | 热处理对豌豆分离蛋白结合典型牛肉香气物质的影响及相互作用机制 |
| Authors | **XU Sijia** (徐思佳)^1, **BI Shuang** (毕爽)^2, **ZHANG Wentao** (张文涛)^1, **PAN Xin** (潘鑫)^1, **LAO Fei** (劳菲)^1, **SHEN Qun** (沈群)^1, **WU Jihong** (吴继红)^1,* — ^1 National Engineering Research Center for Fruit and Vegetable Processing, Key Laboratory of Fruit and Vegetable Processing of Ministry of Agricultural and Rural Affairs, Key Laboratory of Non-thermal Food Processing, College of Food Science & Nutritional Engineering, China Agricultural University, Beijing 100083; ^2 College of Food and Health, Beijing Technology and Business University, Beijing 100048. Corresponding: WU Jihong, wjhcau@hotmail.com (ORCID 0000-0001-7895-8305). First author ORCID 0000-0003-1992-2817 |
| Venue | 食品科学 *Food Science* (China), **2022, Vol. 43, No. 19, pp. 28-35**. Article number 1002-6630(2022)19-0028-08; CLC TS214.9; document code A. **"(in Chinese with English abstract)"**, as the paper's own citation block states |
| DOI | **`10.7506/spkx1002-6630-20220613-128`** — printed at the end of the English abstract and repeated twice in the bilingual citation block on p. 29, exactly as given here. Journal site printed as `http://www.spkx.net.cn` |
| Dates | 收稿日期 (received) **2022-06-13** |
| Funding | 十三五 National Key R&D Plan key special project **2017YFD0401202** — **the same grant number that funds Bi 2022** |
| Protein | **Pea protein isolate made in-house** from *Pisum sativum* L. cv. **'Zhongwan' (中豌) series, 2021 harvest**, Shanxi Dongfangliang Life Technology Co. **Protein content by Kjeldahl: (90.15 ± 1.20) %** — a figure Bi 2022 never reported for its own isolate |
| The four ligands | **2-methylpyrazine** (≥99 %, Sigma-Aldrich Shanghai), **2-methyl-3-mercaptofuran** = 2-methyl-3-furanthiol (≥95 %, Sigma-Aldrich), **(E,E)-2,4-decadienal** (≥89 %, Sigma-Aldrich), **5-(2-hydroxyethyl)-4-methylthiazole** (>98 %, TCI). ANS ≥97 % (Sigma-Aldrich). All chosen as **key aroma compounds of cooked beef** after Song Ze 2019 (ref. [3]) |
| Naming | "结合率" = **binding rate**, %, by headspace depletion (Eq. 1); "n" = number of binding sites; "K" = binding constant, L/mol; "nK" = overall binding affinity, L/mol; "h" = **Hill coefficient** describing cooperativity between a single site and the ligand (h = 1 reduces Eq. 6 to the Klotz equation); "Ksv" = Stern-Volmer quenching constant, L/mol; "Ka" = effective quenching constant, L/mol; "表面疏水性" = surface hydrophobicity; "ζ电位" = zeta potential; "粒径" = particle size |
| Companions on disk | **`bi2022_extraction.md` — the direct sibling**: same laboratory, same grant, same buffer, same in-house pea isolate, same Klotz/SPME family, and its author Bi Shuang is second author here; it is cited as ref. [7]. **`guo2019_extraction.md`** — cited as ref. [16] and its numbers quoted in this paper's discussion. **`anantharamkrishnan2020b_extraction.md`** — cited as ref. [20], and it is the source of `matrix_sites.py`'s hexanal-BLG rate bracket. Also `damodaran1981_extraction.md` (ref. [18]), `wang2015_extraction.md` / `wang2022_extraction.md` families |

## 1. Why it matters

**1. It is the pea-protein experiment at process temperature, from the laboratory whose 37 C work the registry
just shipped.** Wave B26 added the first plant-protein rows to `REVERSIBLE_BINDING` from Bi 2022 — three pea
constants at **37 C** — and its own note draws the boundary: "37 C is an in-mouth temperature, not a process
one: nothing here licenses a pea binding constant at 90 or 140 C." **This paper is the same group's follow-up
and it goes to 120 °C.** Its buffer (0.01 mol/L K2HPO4-KH2PO4, pH 7.6), its isolate (in-house alkaline
extraction / isoelectric precipitation from Chinese peas), its methanol aroma stocks (0.01 mol/L), its SPME
fibre (50/30 µm DVB/CAR/PDMS), its GC-MS (Agilent, DB-WAX 30 m × 0.25 mm × 0.25 µm) and its Klotz treatment are
all Bi 2022's, so a Xu number and a Bi number sit on nearly the same experimental footing. **It even fixes Bi's
worst gap: the protein content is measured and printed, (90.15 ± 1.20) %.**

**2. Its numbers convert to the registry's per-gram form with no molar mass needed — which Bi's do not.**
Bi 2022's Klotz K is per mole of pea protein and the dossier's Flags 1 records that the molar mass "is not
printed anywhere in the paper", so `K_g` could only be reached through the phase-ratio-variation pair. **Xu
prints a bound FRACTION and a protein LOADING**, and those two are exactly what
`K_g = (K_water/K_matrix − 1)/protein_g_per_L` needs. The binding rate is defined as
**1 − (headspace peak area with protein / headspace peak area without protein)**, i.e. its complement *is*
K_matrix/K_water, so **K_g = R/((1 − R) × C)** with no molar mass anywhere in the chain (mine, section 3):

| compound | protein state | binding rate printed | **K_g at 5 g/L solids (mine)** |
|---|---|---:|---:|
| 2-methylpyrazine | native | 17.58 ± 0.58 % | **4.27e-2 L/g** |
| 2-methylpyrazine | 100 °C, 10 min | 25.35 ± 1.26 % | **6.79e-2 L/g** |
| 2-methylpyrazine | 120 °C, 10 min | ~49.6 % (from the printed "+ about 32 %") | **~1.97e-1 L/g** |
| 5-(2-hydroxyethyl)-4-methylthiazole | native | 80.34 ± 7.68 % | **8.17e-1 L/g** |
| 5-(2-hydroxyethyl)-4-methylthiazole | 100 °C, 10 min | 54.41 ± 4.14 % | **2.39e-1 L/g** |
| 2-methyl-3-furanthiol, (E,E)-2,4-decadienal | native and heated | **> 98 %** | **> 9.8 L/g** — and see point 4 |

For scale, the shipped `kg_hexanal_pea` is **2.537e-1 L/g** and `kg_t_2_octenal_pea` (quarantined) is
3.834e-1. **2-Methylpyrazine on native pea is 6x weaker than hexanal; after a 120 °C cook it is within 1.3x of
it.** The thiazole starts 3.2x *above* hexanal and falls below it after 100 °C.

**3. The preheat effect is bidirectional again, and the two directions are in the same table.**
`matrix_sites.py` charges its site pools once at the start of the cook and does not change them with heating.
Guo 2019 showed the change is ligand-dependent on soy; **this paper shows the same on pea, with the split
running along a different axis**: the compound with the **lowest log P** (2-methylpyrazine, 0.785) **gains**
binding as the protein is heated (+7 percentage points at 100 °C, **+ about 32 points at 120 °C**), while the
sterically bulky thiazole **loses** it (−26 points at 100 °C) before partially recovering at longer times. The
paper's mechanism is explicit and testable: **heating unfolds the protein and exposes hydrophobic surface, and
2-methylpyrazine's sites are ON that surface** — the binding rate correlates with heat-treatment temperature at
**r = 1.00**, with heat-treatment time at **r = 0.90** and with particle size at **r = −0.96** (p. 33). For the
thiazole the correlation with surface hydrophobicity is **r = −0.89**, and the reason offered is competition:
heating raises the compound's activity coefficient and drives it into the headspace faster than the protein can
bind it. **So the model's "preheat multiplier" is not a scalar and is not even one-signed on pea.**

**4. The two strongest binders are declared COVALENT by the authors, which quarantines them cleanly.** This is
unusually helpful. Section 2.2 (p. 32) says of 2-methyl-3-furanthiol that it "contains a thiol and can interact
strongly with protein", citing Anantharamkrishnan 2020 for **disulfide and trisulfide adducts** with
beta-lactoglobulin; and of (E,E)-2,4-decadienal that as an unsaturated aldehyde it forms **irreversible covalent
bonds**, **Schiff bases with the protein's amide side chains** and **conjugate (Michael) addition to free
cysteine thiols**, and that **heat treatment promotes these reactions**. It then states the consequence in one
sentence: **"由于共价结合是不可逆的" — because covalent binding is irreversible** — flavour is lost and shelf
life shortened, and only non-covalent interaction preserves flavour through processing. **The >98 % binding of
these two compounds is therefore NOT a reversible binding constant on the paper's own account, and both must
be refused entry to `REVERSIBLE_BINDING`** — the same ruling `kg_t_2_hexenal_dairy` (Meynier) and
`kg_t_2_octenal_pea` (Bi) already carry, but here made by the authors rather than inferred. **(E,E)-2,4-decadienal
is precisely the 2-alkenal hazard class**, and **2-methyl-3-furanthiol is a thiol feeding
`matrix_sites.py`'s `PROT_SS` thiol-to-disulfide exchange channel**, whose pea disulfide density
(0.0257 mmol/g) `protein_matrices.yml` already charges.

**5. It is the corpus's first `dose_added_pre_cook` binding measurement on a plant protein — and that is both
its strength and its central defect.** `parameters_matrix.py` already uses that classification: Brewer 1995's
beef observation is **excluded** from the unsaturation fit because it is "declared HOLD-OUT and reclassified
`dose_added_pre_cook`". **Xu's design is exactly that**: the aroma compounds are pipetted into the vial, the
vial is capped, and *then* it is heated to 80-120 °C. So what the GC sees afterwards is the net of binding,
covalent capture, thermal degradation of the ligand, Maillard consumption of the ligand and — for
2-methylpyrazine, which is **itself a Maillard product**, as the paper says on p. 33
("同时也是美拉德反应的主要产物", "it is also a major product of the Maillard reaction") — **possible
formation of the ligand during the cook.** The control that would separate these is a heated
protein-free blank, and the binding-rate denominator ("香气物质峰面积 without pea protein isolate") is exactly
that, **provided the blank was heated too** — which the Methods do not say (Flags 1). **This is the most
process-realistic design in the batch and the one with the most channels between the dose and the readout.**

**6. On the B26 temperature limit: this paper does NOT lift it, and the distinction is the same as Guo's.**
Every binding rate and every constant here is read out after the vial is **cooled in ice water and shaken at
25 °C for 12 h**. **The measurement temperature is 25 °C.** What has been to 120 °C is the *system*, not the
equilibrium. The correct reading is: **this paper measures what a pea matrix holds onto at eating temperature
after a 120 °C cook in which the aroma was present** — which is arguably the question a plant-meat model
actually asks — and it says **nothing** about the binding constant *at* 120 °C. A `preheat_c` provenance key
(as proposed in `guo2019_extraction.md`) would carry it correctly; putting 120.0 in `temperature_c` would be a
category error.

**7. It supplies the thermodynamics — and they are van 't Hoff quantities, not activation energies.** Table 2
prints **ΔH = −86.68 kJ/mol** (native) and **−60.83 kJ/mol** (100 °C-treated) for 2-methylpyrazine × pea. These
are fitted from the temperature dependence of the **fluorescence-quenching constant Ka** over 298/303/310 K by
Eq. 4, ln Ka = −ΔH/(RT) + ΔS/R. **They are enthalpies of a binding equilibrium. They are NOT activation
energies**, and −86.68 kJ/mol is a number that would look entirely plausible dropped into
`matrix_sites.py`'s `ea_band_kj_mol` field (which carries 15-20 kJ/mol for the aldehyde-amine channel, from
real rate measurements) and would be wrong by a factor of 4-6 with the wrong sign convention besides. **An
activation energy cannot be negative; this ΔH is.** Note also the paper prints the gas constant as
"8.314 kJ/(mol·K)" (p. 30) — a unit slip; the computation used J/(mol·K), which this dossier verifies
(section 3).

**What this paper does NOT give the repository**: any binding measurement above 25 °C (the fluorescence work
reaches 37 C); any molar mass for pea protein, so the Klotz K in L/mol cannot be converted (Flags 5); any rate
constant or activation energy; any adduct identification or mass-spectrometric adduct search (the covalent
claim is by citation, not by measurement here); any binding constant for the thiazole, the furanthiol or the
decadienal (only 2-methylpyrazine gets a Klotz fit); any pH other than 7.6; any protein-free heated blank that
is described as such; any numeric binding rate other than the five quoted; any sensory measurement or odour
threshold.

## 2. Methods as they matter to a model

- **The protein.** *Pisum sativum* L. cv. 'Zhongwan' series, 2021 harvest, Shanxi Dongfangliang. After Cui Leqi
  2020 with modifications: milled and passed a **60-mesh** sieve; ultrapure water at **1:10 (m/V)**; pH raised
  to **9.0 with 2 mol/L NaOH**; stirred **1 h**; **10 000 × g, 15 min**; supernatant taken to **pH 4.5 with
  2 mol/L HCl**; stirred **1 h**; **10 000 × g, 15 min**; supernatant discarded and the precipitate washed
  **twice** with ultrapure water; pH raised to **7.0 with 2 mol/L NaOH**; transferred to a **7 kDa dialysis
  bag** and dialysed **48 h at 4 °C with the water changed every 6 h**; freeze-dried. **Protein content by
  Kjeldahl (90.15 ± 1.20) %.** This is the same preparation family as Bi 2022's, with a dialysis step Bi did not
  have and a protein assay Bi did not report. **No thiol, disulfide, amine or molar-mass determination.**
- **The pot.** **0.01 mol/L K2HPO4-KH2PO4 buffer, pH 7.6** throughout — identical to Bi 2022. A **10 mg/mL**
  pea protein solution was made in it and magnetically stirred **12 h** to dissolve fully. Aroma stocks
  **0.01 mol/L in chromatographic-grade methanol**, diluted with the same phosphate buffer **freshly before
  each experiment** to **40 mg/L** (2-methylpyrazine, 2-methyl-3-furanthiol, (E,E)-2,4-decadienal) and
  **220 mg/L** (5-(2-hydroxyethyl)-4-methylthiazole).
- **Loading — read this carefully, it is HALF Bi 2022's.** Into a **20 mL headspace vial**: **1 mL of the
  10 mg/mL protein solution + 0.5 mL of the diluted aroma + 0.5 mL phosphate buffer = 2 mL**. So the **final
  protein is 5 mg/mL = 5 g/L of isolate**, and at 90.15 % protein, **4.51 g/L of protein (mine)**. Bi 2022's
  headspace assay ran at 10 g/L. **The final aroma concentrations (mine, 1:4 dilution of the stated stocks) are
  10 mg/L = 0.106 mmol/L (2-methylpyrazine), 10 mg/L = 0.088 mmol/L (2-methyl-3-furanthiol), 10 mg/L =
  0.066 mmol/L ((E,E)-2,4-decadienal) and 55 mg/L = 0.384 mmol/L (the thiazole).**
- **THE HEAT TREATMENT — and the aroma is inside the vial for it.** "混合后立即旋紧螺旋盖。然后在各热处理
  条件下处理样品" — *after mixing, the screw cap was immediately tightened; the samples were then treated under
  each heat-treatment condition.* **Temperature arm: 80, 90, 100, 120 °C for 10 min. Time arm: 100 °C for 5, 10,
  15, 20 min. Control: no heat treatment.** Immediately after heating, **the vial was put into ICE WATER to cool**,
  then **shaken at a constant 25 °C for 12 h to reach equilibrium**, and only then sampled. **So the ligand is
  present through the entire thermal excursion, and the binding equilibrium is established afterwards at 25 °C.**
- **Measurement family: HEADSPACE SPME-GC-MS depletion.** The vial was **equilibrated at 50 °C for 10 min with
  shaking**, then the SPME fibre (50/30 µm DVB/CAR/PDMS, 2 cm) was exposed to the headspace for **30 min**, and
  desorbed at **250 °C for 5 min**. GC: helium (99.999 %) at **1.0 mL/min**, **1 mL** of headspace drawn into
  the inlet, **splitless**, inlet 250 °C, oven **40 °C hold 2 min → 10 °C/min → 220 °C hold 2 min**. MS: EI,
  **70 eV**, ion source 230 °C, auxiliary heater 250 °C, quadrupole 150 °C. **Quantitation in SIM.**
  Characteristic ions determined in full scan: **2-methylpyrazine m/z 94, 67, 26, 39, 40**;
  **2-methyl-3-furanthiol m/z 114, 113, 85, 45, 43**; **(E,E)-2,4-decadienal m/z 81, 41, 67, 83, 55**;
  **5-(2-hydroxyethyl)-4-methylthiazole m/z 112, 113, 45, 143, 85**.
  **Note the sampling temperature is 50 °C, not the 25 °C at which the binding equilibrium was set** (Flags 6).
  For the registry's `method` field this is **`headspace_depletion`** — the same class as Bi 2022's Klotz K and
  Andriot 2000's beta-lactoglobulin rows, and firmly on the headspace side of the k2 sec. B.3 boundary.
- **Equation 1, the binding rate** (recovered from the page image; the text layer rendered it as mojibake):
  **结合率 / % = (1 − [peak area of the aroma compound in the sample CONTAINING pea protein isolate] / [peak
  area of the aroma compound in the sample NOT containing pea protein isolate]) × 100**, after Wang Kun &
  Arntfield 2017 (ref. [6]). **This is a two-sample ratio at fixed total dose**, which is why it converts
  directly into the registry's `K_g` form (section 3).
- **Surface hydrophobicity (ANS).** Pea protein at **0.005, 0.010, 0.015, 0.020, 0.025, 0.030 mg/mL** in the
  same buffer, **heat-treated**, then 4 mL + **20 µL of 8 mmol/L ANS**, held **5 min at 25 °C in the dark**;
  excitation **390 nm**, emission **470 nm**, spectrum 430-520 nm, slit 1 nm, resolution 0.5 nm. **Surface
  hydrophobicity = the slope of fluorescence intensity against protein concentration.**
- **Particle size and zeta potential.** After heat treatment, diluted to **1 mg/mL** in the same buffer; Malvern
  ZEN 3700; **relative refractive index 1.330**.
- **Fluorescence quenching.** Pea protein **1 mg/mL**; 2-methylpyrazine at **0, 0.04, 0.08, 0.12, 0.16,
  0.20 mmol/L**; made to **5 mL** with the same buffer. **Control = unheated; heat-treated group = 100 °C for
  10 min.** Then shaken at **100 r/min for 2 h at 25, 30 or 37 °C**, and the fluorescence measured immediately.
  Excitation **290 nm**, emission **330 nm**, spectrum 300-450 nm, slit 2 nm, resolution 0.5 nm. Stern-Volmer
  **F0/F = 1 + Ksv[Q]** (Eq. 2) gives Ksv; the modified form **F0/ΔF = F0/(F0 − F) = 1/(fa·Ka·[Q]) + 1/fa**
  (Eq. 3) gives the effective quenching constant Ka; **ln Ka = −ΔH/(RT) + ΔS/R** (Eq. 4) gives ΔH and ΔS; and
  **ΔG = ΔH − TΔS** (Eq. 5). **The paper prints R as "8.314 kJ/(mol·K)" — a unit slip; the arithmetic used
  J/(mol·K), which this dossier verifies (section 3 item 4).** Note the protein loading here is **1 mg/mL, five
  times lower than the 5 mg/mL of the headspace assay**.
- **The binding-constant experiment (§1.3.6).** In a 20 mL headspace vial: **pea protein at 5 mg/mL**,
  **2-methylpyrazine at 0.06, 0.08, 0.10, 0.12, 0.14 mmol/L**, made to **2 mL** with the same buffer.
  **Control = unheated; heat-treated group = 100 °C for 10 min.** SPME-GC-MS as in §1.3.2. The Scatchard model
  (ref. [13]) with **Eq. 6: 1/v = 1/n + (1/n)·(K[L])^(−h)** and **Eq. 7: v = ([L]t − [L])/Cp**, where v is
  moles of volatile bound per mole of protein, [L]t the total volatile concentration, **Cp the pea protein
  concentration in mol/L**, n the number of sites, K the binding constant in L/mol, [L] the free concentration
  and **h the Hill coefficient describing cooperativity between a single binding site and the ligand; when
  h = 1, Eq. 6 is the Klotz equation, in which nK is customarily used to measure overall binding affinity.**
  **The fitted value of h is never printed** (Flags 5). **Cp requires a molar mass and none is given.**
- **Statistics.** All experiments in **triplicate**; results as mean ± standard deviation; SPSS v.26.0 for
  one-way ANOVA (P < 0.05), **analysis of covariance** and **Spearman correlation**; Origin 2022 for plotting.

## 3. Tables re-typed

All three tables were verified against the rendered page images. Chinese terms are given with an English
translation; the numbers are as printed. Evidence marks: `[M]` measured in this study, `[C]` cited from
elsewhere, `[F]` fitted.

### Table 1 (p. 32). 表1 4 种典型牛肉香气物质的性质 / "Properties of four typical beef aroma substances"

| 香气物质 (aroma substance) | CAS号 (CAS) | 分子式 (formula) | 相对分子质量 (rel. molecular mass) | log P | 感官描述词 (sensory descriptor) |
|---|---|---|---:|---:|---|
| 2-甲基吡嗪 **2-methylpyrazine** | 109-08-0 `[C]` | C5H6N2 `[C]` | 94.11 `[C]` | **0.785** `[C]` | 烧烤及坚果味 — *roasted / grilled and nutty* |
| 2-甲基-3-巯基呋喃 **2-methyl-3-mercaptofuran** (2-methyl-3-furanthiol) | 28588-74-1 `[C]` | C5H6OS `[C]` | 114.17 `[C]` | **1.877** `[C]` | 肉香 — *meaty* |
| (E,E)-2,4-癸二烯醛 **(E,E)-2,4-decadienal** | 25152-84-5 `[C]` | C10H16O `[C]` | 152.23 `[C]` | **2.878** `[C]` | 脂肪味、鸡油味 — *fatty, chicken-fat* |
| 5-羟乙基-4-甲基噻唑 **5-(2-hydroxyethyl)-4-methylthiazole** | 137-00-8 `[C]` | C6H9NOS `[C]` | 143.21 `[C]` | **0.986** `[C]` | 肉味、煮肉、烤肉味 — *meaty, boiled meat, roast meat* |

Footnote exactly as printed: *注：log P表示油水分配系数或疏水系数。* — "log P denotes the oil-water partition
coefficient, or hydrophobicity coefficient."

**Every cell in this table is `[C]`, and its source is a commercial web database.** §1.3.2 states:
*"香气物质自身的属性参数查阅于 https://www.chemicalbook.com/ProductIndex.aspx"* — the aroma compounds' own
property parameters were looked up at chemicalbook.com. **Nothing in Table 1 was measured**, and the log P
column in particular is a database value with no stated method or reference. `parameters_matrix.py` refuses
"any shipped matrix term that is a monotone function of log P" under k4b hold-out guard #4; **this table cannot
be used to build one, and the correlation the paper draws from it (section 2.2) rests on database values, not
measurements** (Flags 3).

### Table 2 (p. 34). 表2 热处理前后2-甲基吡嗪与豌豆分离蛋白相互作用的Stern-Volmer常数和热力学参数 / "Stern-Volmer constants and thermodynamic parameters of interaction between 2-methylpyrazine and pea protein isolate before and after heat treatment"

| 样品 (sample) | T / K | Ksv / (L/mol) | Ka / (L/mol) | ΔH / (kJ/mol) | ΔS / (kJ/(mol·K)) | ΔG / (kJ/mol) |
|---|---:|---|---|---:|---:|---:|
| **对照组** (control, unheated) | 298 | 1 937.62 ± 103.08 ^a `[M]` | 4 220.15 ± 188.49 ^b `[F]` | **−86.68** `[F]` | **−0.22** `[F]` | −20.65 `[F]` |
| | 303 | 1 720.11 ± 54.24 ^b `[M]` | 2 279.02 ± 183.97 ^d `[F]` | | | −19.54 `[F]` |
| | 310 | 1 467.46 ± 34.86 ^d `[M]` | 1 085.4 ± 111.39 ^f `[F]` | | | −17.99 `[F]` |
| **热处理组** (heat-treated, 100 °C / 10 min) | 298 | 1 725.15 ± 105.22 ^b `[M]` | 4 772.41 ± 248.72 ^a `[F]` | **−60.83** `[F]` | **−0.13** `[F]` | −20.77 `[F]` |
| | 303 | 1 682.26 ± 38.56 ^b `[M]` | 2 512.65 ± 38.42 ^c `[F]` | | | −20.10 `[F]` |
| | 310 | 1 520.89 ± 16.17 ^c `[M]` | 1 804.48 ± 45.53 ^e `[F]` | | | −19.16 `[F]` |

Footnote exactly as printed: *注：同列肩标小写字母不同表示差异显著（P＜0.05）。表3同。* — "different lower-case
superscripts **in the same column** indicate a significant difference (P < 0.05); the same applies to Table 3."
(Note that this footnote is correct about its own layout, unlike Guo 2019's.)

Fit qualities stated in the text (p. 33): Stern-Volmer **R^2 > 0.98**; the modified Stern-Volmer form
**R^2 > 0.99**; the van 't Hoff fit **R^2 > 0.93**.

**Ksv falls with rising temperature in both groups, which the paper reads as STATIC quenching** (p. 33) — a
ground-state complex, not a diffusive collision. **ΔH < 0 and ΔS < 0 in both groups, which by Ross &
Subramanian 1981 (ref. [27]) assigns the interaction to van der Waals forces and hydrogen bonds**, and the
paper calls it "焓驱动" — enthalpy-driven.

**Internal check (mine): does the van 't Hoff arithmetic close?** Refitting ln Ka against 1/T over the three
printed Ka values with **R = 8.314 J/(mol·K)**: the control gives **ΔH = −86.68 kJ/mol and ΔS = −0.2216
kJ/(mol·K)** against a printed −86.68 and −0.22 — **exact**; the heat-treated group gives **ΔH = −60.83 and
ΔS = −0.1344** against a printed −60.83 and −0.13 — **exact**. **So the printed gas constant "8.314 kJ/(mol·K)"
on p. 30 is a unit slip and the computation used J/(mol·K) correctly.**

**Internal check (mine): ΔG against −RT ln Ka.** Control: −20.68, −19.48, −18.01 kJ/mol at 298/303/310 K
against printed −20.65, −19.54, −17.99 — **agreement to 0.5 %**. Heat-treated: −20.99, −19.72, −19.32 against
printed −20.77, −20.10, −19.16 — **agreement to 2 %**. Both close; the residual is the rounding of ΔS to two
decimals in the display.

### Table 3 (p. 34). 表3 热处理前后2-甲基吡嗪与豌豆分离蛋白的结合参数 / "Binding parameters of 2-methylpyrazine to pea protein isolate before and after heat treatment"

| 组别 (group) | n | K / (L/mol) | nK / (L/mol) |
|---|---|---|---|
| **对照组** (control, unheated) | 20.46 ± 15.69 ^a `[F]` | 638.63 ± 490.36 ^a `[F]` | **13 065.25 ± 546.63 ^b** `[F]` |
| **热处理组** (heat-treated, 100 °C / 10 min) | 23.78 ± 15.12 ^a `[F]` | 961.97 ± 612.93 ^a `[F]` | **22 876.70 ± 985.48 ^a** `[F]` |

Same footnote as Table 2 (letters within a column). Fig. 5's Klotz plot carries **R^2 = 0.979** (pea protein
isolate) and **R^2 = 0.978** (heat-denatured pea protein isolate) inside the panel.

**Read the letters: n and K are NOT significantly different between the two groups (both ^a); only nK is
(^b against ^a).** This is exactly what the error bars say — **the standard error on n is 77 % of its value and
on K is 77 % and 64 %, while on nK it is 4.2 % and 4.3 % (mine)** — the signature of a double-reciprocal fit in
which the slope is well determined and the intercept is not. **The only defensible number in this table is nK.**

**Internal check (mine): n × K against the printed nK.** Control: 20.46 × 638.63 = **13 066.4** against
13 065.25 ✓ (0.01 %). Heat-treated: 23.78 × 961.97 = **22 875.6** against 22 876.70 ✓ (0.005 %). **The table
closes to five significant figures.**

### Numbers printed in the running text

| quantity | value | where | class |
|---|---|---|---|
| protein content of the isolate | **(90.15 ± 1.20) %** by Kjeldahl | §1.3.1, p. 30 | `[M]` |
| protein stock / final loading | **10 mg/mL** stock; 1 mL + 0.5 mL aroma + 0.5 mL buffer in a 20 mL vial | §1.3.2, p. 30 | `[M]` — **implies 5 g/L isolate, 4.51 g/L protein (mine)** |
| aroma stock | **0.01 mol/L in methanol**, diluted in buffer to **40 mg/L** (2-MP, MFT, decadienal) and **220 mg/L** (thiazole) | §1.3.2, p. 30 | `[M]` |
| **heat-treatment temperature arm** | **80, 90, 100, 120 °C for 10 min**, aroma present, then ice water, then 12 h at 25 °C | §1.3.2, p. 30 | `[M]` |
| **heat-treatment time arm** | **100 °C for 5, 10, 15, 20 min**, otherwise as above | §1.3.2, p. 30 | `[M]` |
| SPME | equilibrate **50 °C, 10 min**; adsorb **30 min**; desorb **250 °C, 5 min** | §1.3.2, p. 30 | `[M]` |
| **binding rate, 2-methyl-3-furanthiol and (E,E)-2,4-decadienal** | **> 98 %**, at room temperature AND after heating; "热处理对结合率的影响不大" (heat treatment has little effect) | Abstract; §2.1.1, p. 31; §3, p. 34 | **`[M]`** |
| **binding rate, 5-(2-hydroxyethyl)-4-methylthiazole, room temperature** | **(80.34 ± 7.68) %** | Abstract; §2.1.1, p. 31 | **`[M]`** |
| **binding rate, 2-methylpyrazine, room temperature** | **(17.58 ± 0.58) %** | Abstract; §2.1.1, p. 31 | **`[M]`** |
| **binding rate after 100 °C / 10 min, thiazole and 2-methylpyrazine** | **(54.41 ± 4.14) %** and **(25.35 ± 1.26) %** | Abstract; §2.1.2, p. 31 | **`[M]`** |
| 2-methylpyrazine at 100 °C vs room temperature | "增加了约 7 %" — up by **about 7 percentage points** | §2.1.1, p. 31 | `[M]` (see the note below) |
| thiazole at 100 °C vs room temperature | "降低了约 26 %" — down by **about 26 percentage points** | §2.1.1, p. 31 | `[M]` |
| **2-methylpyrazine at 120 °C / 10 min vs room temperature** | **"大幅提高了约 32 %"** — up by **about 32 percentage points** | §2.1.1, p. 31 | **`[M]` — the highest-temperature binding datum in the paper** |
| thiazole at 120 °C / 10 min | "略低于常温下的结合率，但是高于 80～100 °C 的结合率" — slightly below the room-temperature rate but above the 80-100 °C rates | §2.1.1, p. 31 | `[M]` (no value) |
| 2-methylpyrazine after 15 min at 100 °C | "相比常温时的结合率增加了约 12 %" — up **about 12 points** on room temperature | §2.1.2, p. 31 | `[M]` |
| 2-methylpyrazine, 20 min vs 15 min | **no significant difference (P > 0.05)** — "说明结合已经达到平衡", binding has reached equilibrium | §2.1.2, p. 31 | `[M]` |
| thiazole after 5 min at 100 °C | "相比常温时下降了约 44 %" — down **about 44 points** | §2.1.2, p. 31 | `[M]` |
| thiazole after 20 min at 100 °C | "上升到略高于常温" — risen to slightly above the room-temperature value | §2.1.2, p. 31 | `[M]` (no value) |
| 120 °C verdict | "表明 120 °C 处理能够促进豌豆分离蛋白结合牛肉香气物质" — 120 °C treatment **promotes** the binding of beef aroma substances by pea protein isolate | §2.1.1, p. 31 | `[M]`/interpretation |
| **correlation, 2-methylpyrazine binding rate vs heat-treatment temperature** | **r = 1.00** | §2.3.2, p. 33 | **`[F]`** (Spearman) |
| correlation, 2-methylpyrazine binding rate vs heat-treatment time | **r = 0.90** | §2.3.2, p. 33 | `[F]` |
| correlation, 2-methylpyrazine binding rate vs particle size | **r = −0.96** | §2.3.2, p. 33 | `[F]` |
| **correlation, thiazole binding rate vs surface hydrophobicity** | **r = −0.89** | §2.3.1, p. 33 | `[F]` |
| surface hydrophobicity, temperature arm | rises to a **maximum at 100 °C**; **falls slightly at 120 °C** — "表明疏水残基不再暴露，倾向于聚集" (hydrophobic residues no longer exposed, tending to aggregate) | §2.3.1, p. 32, Fig. 3A | `[M]` — **magnitudes are figure-only** |
| surface hydrophobicity, time arm | falls slightly with prolonged heating at 100 °C | §2.3.1, p. 32, Fig. 3B | `[M]` — figure-only |
| particle size | **falls significantly (P < 0.05)** on heating and keeps falling as the temperature rises — "表明热处理诱导了蛋白四级结构的解离" (heat treatment induces dissociation of the quaternary structure); **no significant change with heating TIME** | §2.3.2, p. 32, Fig. 4A-B | `[M]` — figure-only |
| zeta potential | rises slightly with temperature and time but **stays near −30 mV**, so "体系相对稳定，蛋白质不会大量聚集" (the system is relatively stable and the protein does not aggregate substantially) | §2.3.2, p. 32, Fig. 4C-D | `[M]` — figure-only |
| the stability criterion cited | a protein is stable in solution when \|ζ\| is not below **30 mV** | §2.3.2, p. 32 | `[C]` (Hartmann & Palzer 2011) |
| 11S behaviour at 80 °C | not fully denatured but "结构变得更加紧密" (structure becomes more compact) to adapt to the rising temperature | §2.3.2, p. 32 | `[C]` (Sorgentini 1995) |
| why no aggregation | at low ion concentration and low ionic strength free subunits struggle to form aggregates; denatured 7S and convicilin reduce 11S aggregation by steric hindrance and electrostatic repulsion, both non-covalent | §2.3.2, p. 32 | `[C]` (Guo 2019; Mession 2013) |
| **the covalent statement, 2-methyl-3-furanthiol** | contains a thiol; beta-lactoglobulin forms **disulfide and trisulfide adducts** with propanethiol, furfuryl mercaptan and thiophenol **by covalent bonds** | §2.2, p. 32 | **`[C]` (Anantharamkrishnan 2020, ref. [20])** |
| **the covalent statement, (E,E)-2,4-decadienal** | an unsaturated aldehyde; its aldehyde group interacts strongly with protein functional groups forming **irreversible covalent bonds**; alkenals form **Schiff bases with amide side chains** and undergo **conjugate addition with free cysteine thiols**; **heat treatment promotes these reactions** | §2.2, p. 32 | **`[C]` (Wang Juan 2018 ref. [21]; Anantharamkrishnan 2020 ref. [20])** |
| **the consequence the authors draw** | **"由于共价结合是不可逆的"** — because covalent binding is irreversible, it causes aroma loss and shortened shelf life; **only non-covalent interaction preserves flavour through processing** | §2.2, p. 32 | **interpretation — the quarantine ruling, in the authors' own words** |
| why the pyrazine and the thiazole bind weakly at room temperature | both contain heterocycles, are **compact and rigid**, and their **large steric hindrance** limits hydrophobic interaction with the protein; the flexible-chain (E,E)-2,4-decadienal can change its own conformation to reach sites | §2.2, p. 32 | interpretation |
| the log P correlation claimed | binding rate is **positively correlated with log P**; the lowest-log P compound (2-methylpyrazine) binds least | §2.2, p. 32 | interpretation — **on database log P values (Flags 3)** |
| 2-methylpyrazine's identity | "吡嗪类化合物是对炖煮牛肉香气有重要贡献的杂环类化合物，**同时也是美拉德反应的主要产物**" — pyrazines contribute importantly to stewed-beef aroma and **are also major Maillard reaction products** | §2.4, p. 33 | **`[C]`/interpretation — the reason the 120 °C arm needs a synthesis control (Flags 1)** |
| Guo 2019's numbers, as quoted here | soy isolate: high-affinity primary sites **n = 1-2, K = 14 000-20 000 L/mol**; low-affinity secondary sites **n = 1-11, K = 100-1 500 L/mol** | §2.4.2, p. 34 | `[C]` — **and the quotation is loose (Flags 8)** |
| the conclusion drawn from that comparison | 2-methylpyrazine on pea has **more sites and a smaller constant** than the soy values, so its sites are inferred to lie **on the hydrophobic SURFACE** rather than in an internal cavity | §2.4.2, p. 34 | interpretation |

**A translation note that changes the numbers.** The Chinese "增加了约 7 %" / "降低了约 26 %" read literally as
"increased by about 7 %" / "decreased by about 26 %", but they are **percentage-POINT changes**, not relative
ones: 17.58 + 7 = 24.58 against the printed 25.35, and 80.34 − 26 = 54.34 against the printed 54.41 **(mine —
both reconstructions land within 0.8 of the independently printed values, so the percentage-point reading is
confirmed)**. The same reading gives **2-methylpyrazine at 120 °C ≈ 17.58 + 32 = 49.6 %** and **at 15 min ≈
29.6 %**, and **the thiazole at 5 min ≈ 36.3 %**. These three are reconstructions, not printed values.

**Figure-only quantities.** All twenty bars of Fig. 1 (four compounds × five heat-treatment temperatures) and
all twenty of Fig. 2 (four compounds × five times), except the five values quoted above; all ten bars of Fig. 3
(surface hydrophobicity, y-axis in units of ×10^5); all twenty of Fig. 4 (particle size in nm, y-axis 0-300;
zeta potential in mV, y-axis −40 to 0); and every point of the two Klotz lines in Fig. 5. Per house rule none is
typed here.

### Arithmetic on the printed constants (all mine)

**1. The binding rates converted to the registry's per-gram form — and this needs NO molar mass (mine).** The
registry stores `K_g = (K_water/K_matrix − 1) / protein_g_per_L`, where K is an air/matrix partition
coefficient. Eq. 1 defines the binding rate as **R = 1 − A_with / A_without**, where A is a headspace peak area
at fixed total dose, so **A_without/A_with = K_aw,water / K_aw,matrix = 1/(1 − R)** and therefore

> **K_g = (1/(1 − R) − 1) / C = R / ((1 − R) · C)**

with C the protein loading. Two loadings are defensible: **5 g/L** of isolate (what was weighed) or **4.51 g/L**
of protein (5 g/L × 90.15 %, the basis the registry's per-gram convention actually means):

| compound | protein state | R printed | **K_g at 5 g/L isolate** | **K_g at 4.51 g/L protein** |
|---|---|---:|---:|---:|
| 2-methylpyrazine | native | 0.1758 | **4.27e-2 L/g** | **4.73e-2 L/g** |
| 2-methylpyrazine | 100 °C, 10 min | 0.2535 | **6.79e-2 L/g** | **7.53e-2 L/g** |
| 2-methylpyrazine | 120 °C, 10 min | ~0.496 (reconstructed) | **~1.97e-1 L/g** | **~2.18e-1 L/g** |
| 5-(2-hydroxyethyl)-4-methylthiazole | native | 0.8034 | **8.17e-1 L/g** | **9.07e-1 L/g** |
| 5-(2-hydroxyethyl)-4-methylthiazole | 100 °C, 10 min | 0.5441 | **2.39e-1 L/g** | **2.65e-1 L/g** |
| 5-(2-hydroxyethyl)-4-methylthiazole | 100 °C, 5 min | ~0.363 (reconstructed) | **~1.14e-1 L/g** | **~1.27e-1 L/g** |
| 2-methyl-3-furanthiol and (E,E)-2,4-decadienal | native and all heated states | > 0.98 | **> 9.8 L/g** | **> 10.9 L/g** |

**Against the shipped pea rows** (`kg_hexanal_pea` 2.537e-1 L/g, `kg_z_2_penten_1_ol_pea` 4.14e-2,
`kg_t_2_octenal_pea` 3.834e-1, all at 37 C from Bi 2022's PRV pair): **native 2-methylpyrazine (4.27e-2)
lands almost exactly on the alcohol row; after 120 °C it lands between hexanal and the quarantined alkenal.**
The thiazole's native value (8.17e-1) is **larger than any row in `REVERSIBLE_BINDING`** and the two covalent
compounds' bound (>9.8 L/g) is **25x larger than the largest** — which is itself a signal that they are not
measuring the same thing (Flags 2).

**2. The preheat factor on pea, as a ratio (mine).**

| compound | 100 °C / native | 120 °C / native |
|---|---:|---:|
| 2-methylpyrazine, on the binding rate | **1.44x** | **~2.82x** (reconstructed) |
| 2-methylpyrazine, on **K_g** | **1.59x** | **~4.61x** (reconstructed) |
| thiazole, on the binding rate | **0.68x** | slightly below 1 (no value printed) |
| thiazole, on **K_g** | **0.29x** | — |
| 2-methyl-3-furanthiol, (E,E)-2,4-decadienal | **~1.00x** (both >98 % before and after) | ~1.00x |

**The spread at 100 °C runs from 0.29x to 1.59x on K_g across three compounds in one experiment (mine)** — a
factor of 5.5 — and adding the 120 °C pyrazine widens it to **16x**. **This is the pea replication of Guo 2019's
soy finding, on a different compound class, and it says the same thing: a scalar preheat multiplier on a site
pool is refuted.** Note also that the direction here is the *opposite* of Crowther 1980's dry-soy autoclave
result (0.51-0.66x for every compound) — **three papers, three matrices, three signs.**

**3. The Klotz constants as a ratio (mine).** From Table 3: **nK rises 22 876.70 / 13 065.25 = 1.751x** after
100 °C / 10 min; K rises 1.506x and n rises 1.162x, but neither of those two is significant on the paper's own
letters. **Compare the same treatment's effect on the binding rate, 1.44x, and on K_g, 1.59x (item 2).** The
three routes agree to within 22 % on the size of the 100 °C effect — **1.44x (headspace depletion at 0.106 mM),
1.59x (the same, as K_g) and 1.75x (Klotz nK over 0.06-0.14 mM)** — which is a genuinely reassuring internal
cross-check and the strongest thing in this paper.

**4. The van 't Hoff arithmetic verified, and the gas-constant unit slip resolved (mine).** Refitting the three
printed Ka values against 1/T with R = 8.314 **J**/(mol·K) reproduces the printed ΔH and ΔS **exactly** for both
groups (section 3, Table 2). With R in kJ/(mol·K) as the paper prints it, the answers would be 1000x too large.
**The text has a unit typo; the numbers are right.**

**5. The quenching constants and the headspace constants are different objects, and here they differ by ~14x
(mine).** At 310 K the effective quenching constant Ka is **1 085.4 (native)** and **1 804.48 L/mol
(heat-treated)**; the Klotz nK from headspace depletion at 25 °C is **13 065.25** and **22 876.70 L/mol**.
Taking nK against Ka, the ratio is **12.0x (native)** and **12.7x (heat-treated)**. Even against K alone
(638.63 and 961.97) the two disagree by 1.7x and 1.9x in the *other* direction. **Bi 2022's dossier records
exactly this phenomenon on the same protein in the same laboratory — a 5-56x gap between its headspace Klotz K
and its fluorescence Ka — and here it reappears at 12-13x.** The two determinations also run at **5 mg/mL
(headspace) against 1 mg/mL (fluorescence)**, a 5x loading difference on top of the method difference. **The
registry's refusal to cross the `method` boundary is vindicated a second time on this laboratory's own data.**
**Note, however, that the two agree on the DIRECTION and roughly on the SIZE of the preheat effect**: Ka rises
1.66x at 310 K and nK rises 1.75x. That agreement, across two methods with a 12x scale gap, is the best
evidence in the paper that the preheat effect is real.

**6. Ksv falls with temperature in both groups, but by different amounts (mine).** Control: 1 937.62 → 1 467.46
from 298 to 310 K, a **1.32x fall**. Heat-treated: 1 725.15 → 1 520.89, a **1.13x fall**. Both are the static-
quenching signature the paper claims. **The heat-treated protein's quenching is less temperature-sensitive**,
which is consistent with its smaller |ΔH| (60.83 against 86.68 kJ/mol) and is the one place where the Table 2
thermodynamics say something the binding data do not.

**7. The final ligand concentrations, which the paper gives only as mass concentrations of the intermediate
stock (mine).** 0.5 mL of the diluted stock goes into 2.0 mL, a 1:4 dilution: **2-methylpyrazine 10 mg/L =
0.106 mmol/L; 2-methyl-3-furanthiol 10 mg/L = 0.088 mmol/L; (E,E)-2,4-decadienal 10 mg/L = 0.066 mmol/L;
5-(2-hydroxyethyl)-4-methylthiazole 55 mg/L = 0.384 mmol/L.** Note the **thiazole is dosed 3.6x higher in
molar terms than the pyrazine**, so their binding rates are not measured at a common occupancy — the compound
with the higher dose also has the higher bound fraction, which is the opposite of what saturation would predict
and therefore does not explain the ordering, but it does mean the four rates are not strictly comparable
(Flags 4). The Klotz ladder (0.06-0.14 mmol/L) brackets the 0.106 mmol/L used in the binding-rate run, so those
two experiments *are* on a common concentration basis for 2-methylpyrazine.

**8. Methanol is present, and its fraction can be bounded (mine).** The aroma stocks are 0.01 mol/L in
methanol and are diluted with buffer before use. At 40 mg/L of a ~100 g/mol compound the dilution from
0.01 mol/L (≈ 1 000 mg/L) is about **25-fold**, so the diluted aroma solution is roughly **4 % methanol**, and
0.5 mL of it in 2 mL makes the assay about **1 % v/v methanol**. **This is an estimate from the printed
concentrations, not a printed value**, and it is the same order as the ~1.25 % Bi 2022 carries. Methanol
competes for hydrophobic sites and is never mentioned again in either paper.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** **`2_methyl_3_furanthiol` is keyed** (id
`2_methyl_3_furanthiol`, display "2-Methyl-3-furanthiol (MFT)"). **`e_e_2_4_decadienal` is keyed** (and carries
`seen_in: []` — this paper and `zhou2002_extraction.md` are both candidates to populate it).
**`2_methylpyrazine` is ABSENT**, although `2_3_dimethylpyrazine`, `2_5_dimethylpyrazine`,
`2_6_dimethylpyrazine` and `2_ethyl_3_5_dimethylpyrazine` are all present — so the monomethyl parent of a
family the registry already carries four members of is the one missing key. **`5-(2-hydroxyethyl)-4-methylthiazole`
is ABSENT** (`2_hexyl_4_methylthiazole` and `2_pentyl_4_methylthiazole` are present, but not this one).
`COMPOUND_STRUCTURE` in `parameters_matrix.py` carries no pyrazine, no thiazole and no furanthiol, so shipping
either binding row means adding a structural class as well as a key. `MATRIX_LOADING` has no pea entry, but
`pea_protein_1pct` already exists as a `medium` string on the three Wave B26 rows — **note that string says
1 %, i.e. 10 g/L, and this paper is at 5 g/L, so a Xu row needs a distinct medium string.**

Every row below shares: **in-house alkaline-extracted, isoelectric-precipitated, dialysed, freeze-dried pea
protein isolate at (90.15 ± 1.20) % protein; 5 g/L of isolate (4.51 g/L protein) in 0.01 mol/L K2HPO4-KH2PO4,
pH 7.6, ~1 % v/v methanol; 2 mL in a sealed 20 mL headspace vial; aroma present during the heat treatment;
cooled on ice; equilibrated 12 h at 25 C; SPME 50 C / 10 min equilibration + 30 min adsorption; GC-MS in SIM;
triplicate.**

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **binding rate, 2-methylpyrazine, NATIVE pea** | **17.58 ± 0.58** | % | 25 C readout, no heat treatment, 0.106 mM | Abstract; §2.1.1, p. 31 | **`measured_ratio`** |
| **binding rate, 2-methylpyrazine, after 100 °C / 10 min** | **25.35 ± 1.26** | % | as above, aroma present during the cook | Abstract; §2.1.2, p. 31 | **`measured_ratio`** |
| binding rate, 2-methylpyrazine, after 120 °C / 10 min | **~49.6** | % | reconstructed from the printed "+ about 32 percentage points" | §2.1.1, p. 31 (mine) | `derived_assumption` — **reconstruction, not printed** |
| **binding rate, 5-(2-hydroxyethyl)-4-methylthiazole, NATIVE pea** | **80.34 ± 7.68** | % | 25 C readout, no heat treatment, 0.384 mM | Abstract; §2.1.1, p. 31 | **`measured_ratio`** |
| **binding rate, thiazole, after 100 °C / 10 min** | **54.41 ± 4.14** | % | as above | Abstract; §2.1.2, p. 31 | **`measured_ratio`** |
| binding rate, 2-methyl-3-furanthiol and (E,E)-2,4-decadienal | **> 98** | % | native and every heated state | Abstract; §2.1.1, p. 31 | **`measured_bound`** — and **NOT a reversible constant (Flags 2)** |
| **per-gram binding constant, 2-methylpyrazine × native pea** | **4.27e-2** (5 g/L isolate) / **4.73e-2** (4.51 g/L protein) | L/g | 25 C readout, pH 7.6 | R/((1−R)·C) from §2.1.1 (mine) | **`derived_assumption`** — the registry's own K_g form, **no molar mass required** |
| **per-gram binding constant, 2-methylpyrazine × pea after 100 °C / 10 min** | **6.79e-2 / 7.53e-2** | L/g | as above | (mine) | **`derived_assumption`** |
| per-gram binding constant, 2-methylpyrazine × pea after 120 °C / 10 min | **~1.97e-1 / ~2.18e-1** | L/g | as above, from a reconstructed R | (mine) | `derived_assumption` — **doubly derived; treat as indicative** |
| **per-gram binding constant, thiazole × native pea** | **8.17e-1 / 9.07e-1** | L/g | as above | (mine) | `derived_assumption` — **larger than any shipped row** |
| per-gram binding constant, thiazole × pea after 100 °C | **2.39e-1 / 2.65e-1** | L/g | as above | (mine) | `derived_assumption` |
| per-gram bound, furanthiol and decadienal | **> 9.8 / > 10.9** | L/g | as above | (mine) | **`measured_bound` — QUARANTINE, covalent by the authors' own account** |
| **Klotz nK, 2-methylpyrazine × NATIVE pea** | **13 065.25 ± 546.63** | L/mol (per mole of pea protein; **molar mass NOT stated**) | 25 C readout, pH 7.6, 5 g/L, 0.06-0.14 mM | Table 3, p. 34 | **`binding_constant`** |
| **Klotz nK, 2-methylpyrazine × pea after 100 °C / 10 min** | **22 876.70 ± 985.48** | L/mol, as above | as above | Table 3, p. 34 | **`binding_constant`** |
| Klotz K, native / heat-treated | **638.63 ± 490.36 / 961.97 ± 612.93** | L/mol | as above | Table 3, p. 34 | `binding_constant` — **NOT significantly different (both marked ^a); SE is 64-77 % of the value** |
| Klotz n, native / heat-treated | **20.46 ± 15.69 / 23.78 ± 15.12** | mol per mol protein | as above | Table 3, p. 34 | `binding_constant` — **NOT significantly different; SE is 64-77 %** |
| **PREHEAT FACTOR on pea, 2-methylpyrazine, 100 °C / 10 min** | **1.44x** (binding rate), **1.59x** (K_g), **1.75x** (Klotz nK), **1.66x** (Ka at 310 K) | × | four independent routes | §2.1, Tables 2 and 3 (mine) | **`within_study_ratio`** — **the four routes agree within 22 %** |
| **PREHEAT FACTOR on pea, thiazole, 100 °C / 10 min** | **0.68x** (binding rate), **0.29x** (K_g) | × | as above | §2.1.1 (mine) | **`within_study_ratio`** — **the OPPOSITE sign to the pyrazine** |
| spread of the preheat factor across compounds at 100 °C | **0.29x to 1.59x on K_g** (a factor of 5.5); with 120 °C included, **16x** | × | one experiment, one protein, one buffer | (mine) | **`within_study_ratio`** — refutes a scalar preheat term |
| preheat factor, furanthiol and decadienal | **~1.00x** — heat treatment has little effect because both are already >98 % bound | × | as above | §2.1.1, p. 31 | `measured_bound` — **a ceiling effect, not a null** |
| Stern-Volmer constant Ksv, native pea, 298 / 303 / 310 K | **1 937.62 ± 103.08 / 1 720.11 ± 54.24 / 1 467.46 ± 34.86** | L/mol | 1 mg/mL pea, 2 h at each T, fluorescence | Table 2, p. 34 | `binding_constant` — **fluorescence-derived, not headspace** |
| Stern-Volmer constant Ksv, 100 °C-treated pea, 298 / 303 / 310 K | **1 725.15 ± 105.22 / 1 682.26 ± 38.56 / 1 520.89 ± 16.17** | L/mol | as above | Table 2, p. 34 | `binding_constant` — fluorescence |
| effective quenching constant Ka, native, 298 / 303 / 310 K | **4 220.15 ± 188.49 / 2 279.02 ± 183.97 / 1 085.4 ± 111.39** | L/mol | as above | Table 2, p. 34 | `binding_constant` — fluorescence |
| effective quenching constant Ka, 100 °C-treated, 298 / 303 / 310 K | **4 772.41 ± 248.72 / 2 512.65 ± 38.42 / 1 804.48 ± 45.53** | L/mol | as above | Table 2, p. 34 | `binding_constant` — fluorescence |
| **ΔH, ΔS, ΔG (310 K), native pea × 2-methylpyrazine** | **−86.68 / −0.22 / −17.99** | kJ/mol, kJ/(mol·K), kJ/mol | van 't Hoff over 298-310 K on Ka, 1 mg/mL | Table 2, p. 34 | **`binding_constant`** — **a THERMODYNAMIC enthalpy, NOT an activation energy (Flags 9)** |
| ΔH, ΔS, ΔG (310 K), 100 °C-treated pea × 2-methylpyrazine | **−60.83 / −0.13 / −19.16** | as above | as above | Table 2, p. 34 | `binding_constant` — thermodynamic |
| quenching mechanism | **STATIC** in both groups (Ksv falls with T) | — | fluorescence | §2.4.1, p. 33 | `structural_gate` |
| interaction forces assigned | **van der Waals + hydrogen bonds** (ΔH < 0, ΔS < 0), **enthalpy-driven** | — | Ross & Subramanian 1981 | §2.4.1, p. 33 | `derived_assumption` (an assignment rule, not a measurement) |
| method disagreement, headspace nK vs fluorescence Ka at 310 K | **12.0x (native), 12.7x (heat-treated)** | × | same lab, same protein, same buffer | (mine) | **`within_study_ratio`** — the second corroboration of the `method` boundary from this laboratory |
| **the covalent quarantine, in the authors' words** | 2-methyl-3-furanthiol forms **disulfide/trisulfide adducts**; (E,E)-2,4-decadienal forms **Schiff bases and Michael adducts**, **irreversibly**, and **heat promotes both** | — | §2.2, p. 32 | **`structural_gate`** — `[C]`, cited from Anantharamkrishnan 2020 and Wang Juan 2018 |
| protein content of the isolate | **90.15 ± 1.20** | % | Kjeldahl | §1.3.1, p. 30 | `level_only` — **the figure Bi 2022 lacks** |
| Spearman correlations, 2-methylpyrazine | **r = 1.00** (temperature), **0.90** (time), **−0.96** (particle size) | — | across the heat-treatment series | §2.3.2, p. 33 | `derived_assumption` — **on n = 5 points; see Flags 7** |
| Spearman correlation, thiazole vs surface hydrophobicity | **r = −0.89** | — | as above | §2.3.1, p. 33 | `derived_assumption` |
| log P values | **0.785 / 1.877 / 2.878 / 0.986** | — | 2-MP / MFT / decadienal / thiazole | Table 1, p. 32 | **`level_only` — `[C]` from chemicalbook.com; NOT a licence for a log-P term (Flags 3)** |
| molecular masses | **94.11 / 114.17 / 152.23 / 143.21** | g/mol | as above | Table 1, p. 32 | `level_only` — `[C]` |
| CAS numbers | **109-08-0 / 28588-74-1 / 25152-84-5 / 137-00-8** | — | as above | Table 1, p. 32 | `level_only` — `[C]` |
| surface hydrophobicity, particle size, zeta potential | direction only | — | 80-120 °C, 5-20 min | Figs. 3 and 4 | **`figure_only`** |

### Can these be put on the same basis as the shipped binding constants, i.e. converted to K_g in L/g?

**Yes for the binding rates — and this paper is the cleanest conversion in the batch, because the bound
fraction and the loading are both printed and no molar mass is needed. No for the Klotz constants, for the same
reason Bi 2022's could not be converted.**

- **The binding-rate route is exact in form.** Eq. 1's ratio of headspace peak areas with and without protein
  *is* K_matrix/K_water, so `K_g = R/((1 − R)·C)` uses only printed quantities. **This is the same construction
  Amendment 4 blesses for Meynier and Leksrisompong — a within-run ratio in which an absolute headspace
  calibration offset cancels — reached here in one step rather than through a phase-ratio-variation pair.**
  `method` = `headspace_depletion`; `ph_of_measurement` = 7.6; `temperature_c` = **25.0** (the readout, not the
  cook); `medium` = a new pea string at **5 g/L**, distinct from Wave B26's `pea_protein_1pct`.
- **What to use for C.** The registry's per-gram constants mean *per gram of protein*. Bi 2022's shipped rows
  used 10 g/L of a preparation **whose protein content is unknown**, so they are implicitly per gram of isolate.
  **Xu prints 90.15 % protein, so a Xu row can be either.** Ship the 4.51 g/L (per-gram-of-protein) form and
  record the alternative, or ship the 5 g/L form to match Bi's implicit convention — **but do not mix them
  across rows without saying so.** The two differ by 11 %.
- **The Klotz constants cannot be converted.** Eq. 7 divides by **Cp, "豌豆蛋白浓度/(mol/L)"** — the pea protein
  concentration in mol/L — and **no molar mass is printed anywhere in the paper**, exactly as in Bi 2022
  (whose Flags 1 records the same). n = 20.46 sites per mole is itself a clue that a large molar mass was used
  (11S hexamer, ~360 kDa, would give ~14 g/L·mol^-1... but this is speculation and no number is printed).
  **The Klotz K and nK travel only as within-study ratios.** Note that the ratio (1.75x) agrees with the
  molar-mass-free K_g ratio (1.59x) to 10 %, which is the useful cross-check.
- **A further, unusual caution about the Hill coefficient.** Eq. 6 is written with an exponent h and the text
  says "当 h = 1 时公式(6)为 Klotz 方程" — when h = 1, Eq. 6 is the Klotz equation. **The fitted h is never
  printed**, and the table is captioned as Klotz parameters, so h = 1 is presumably assumed — but if it was
  fitted and is not 1, then K's units are not L/mol and nK is not an affinity. **This must be resolved before
  any Klotz number is used** (Flags 5).
- **On the temperature limit.** Every constant here is at a **25 °C readout**. **This paper does not license a
  pea binding constant at 90 or 140 C either.** It licenses something narrower and, for a plant-meat model,
  possibly more useful: **what a pea matrix retains at eating temperature after a real 10-minute cook at up to
  120 °C with the aroma present.** Carry the cook temperature in a `preheat_c` provenance key, never in
  `temperature_c`. And note the design difference from Guo 2019, which preheated the *protein alone* and then
  dosed: **Xu's aroma is in the pot for the cook, which is `dose_added_pre_cook` and is the classification
  `parameters_matrix.py` used to EXCLUDE Brewer 1995 from the unsaturation fit.** The same exclusion logic
  applies here to the two covalent compounds and, more weakly, to 2-methylpyrazine because it is itself a
  Maillard product (Flags 1).

**What goes to `matrix_sites.py`, and it is a structural claim rather than a number.** The two >98 %-bound
compounds map onto two channels that module already models: **(E,E)-2,4-decadienal onto
`unsaturated_aldehyde_amine`** (its `k2_bracket` 5.3-7.9e-5 M^-1 s^-1 at 20 C, from trans-2-hexenal, and note
`BINDING_OF_SPECIES` already maps `DECADIENAL` to that class) and **2-methyl-3-furanthiol onto the
thiol-to-disulfide exchange that charges `PROT_SS`** from `protein_matrices.yml`'s pea disulfide density
(0.0257 mmol/g). **This paper does not supply a rate for either** — no rate constant, no time-resolved adduct
measurement, no activation energy appears anywhere in it. What it supplies is a **measured upper bound on the
combined channel: >98 % of both compounds is removed from the headspace by 4.51 g/L of pea protein at 25 C
without any heating at all**, which is a much larger effect than the pseudo-first-order factor those brackets
would produce over a short thermal programme. **That discrepancy is worth investigating and is not resolved
here** — it may be reversible partition rather than adduction, since a >98 % headspace depletion for a log P
2.878 dienal is also what plain hydrophobic partition into a protein phase would give (Flags 2).

## 5. Flags

1. **THE AROMA IS IN THE VIAL FOR THE COOK, AND ONE OF THE FOUR COMPOUNDS IS A MAILLARD PRODUCT.** The samples
   are mixed, capped and then held at 80-120 °C for 10 min. This is `dose_added_pre_cook` — the classification
   `parameters_matrix.py` uses to EXCLUDE Brewer 1995 from the unsaturation fit. Four processes then compete to
   change the headspace peak area, and the binding rate cannot separate them: **(a) reversible binding, (b)
   covalent capture — which the paper itself says heat promotes, (c) thermal degradation or evaporative loss of
   the ligand, and (d) FORMATION of the ligand.** (d) is not hypothetical: the paper states on p. 33 that
   pyrazines "同时也是美拉德反应的主要产物" — are also major Maillard reaction products — and a pea protein
   solution at 120 °C for 10 min is a Maillard system with amine donors present. **2-Methylpyrazine is the
   compound whose binding "rises" most with heat.** The control that would settle it is the binding-rate
   denominator, "the peak area in the sample NOT containing pea protein isolate" — **but the Methods never say
   whether that protein-free blank was heat-treated alongside the sample.** If it was, (c) largely cancels and
   (d) is partly controlled (though a protein-free blank cannot form pyrazines, since it has no amine donor,
   which would make the cancellation one-sided and inflate the apparent binding). **If it was not, the entire
   temperature series is confounded with simple thermal loss.** This is the first thing to ask the authors.
2. **The two >98 % compounds are covalent by the authors' own account and must be refused as reversible
   constants.** §2.2 (p. 32) states that 2-methyl-3-furanthiol forms disulfide and trisulfide adducts, that
   (E,E)-2,4-decadienal forms irreversible covalent bonds, Schiff bases with amide side chains and Michael
   adducts with free cysteine thiols, that **heat treatment promotes these reactions**, and — in the sentence
   that decides it — **"由于共价结合是不可逆的", because covalent binding is irreversible**. This is the
   `kg_t_2_hexenal_dairy` / `kg_t_2_octenal_pea` precedent stated by the authors rather than inferred by the
   registry. **Neither compound may enter `REVERSIBLE_BINDING`.** Note two further problems with them: the
   binding rates are **at a ceiling** (>98 % before and after heating), so the preheat effect is unmeasurable
   for them, and the >98 % is quoted as an inequality with no error and no value, so **no number exists to
   ship even if the class were right**. Note also that the covalent claim is made **entirely by citation**
   (Anantharamkrishnan 2020, Wang Juan 2018) — **no adduct was identified, no mass spectrum was searched, and
   no reversibility control was run in this paper.** So the compounds are quarantined on cited chemistry, and
   the alternative reading — that a log P 2.878 dienal simply partitions into a protein phase — is untested.
3. **The log P correlation is built on a commercial web database, and the registry refuses log-P terms.** Every
   value in Table 1, including all four log P values, was "looked up at chemicalbook.com" (§1.3.2). No method,
   no primary reference and no uncertainty is given for any of them. The paper's central structural claim —
   "豌豆分离蛋白与4 种牛肉香气物质的结合率与 log P 值呈正相关", the binding rate is positively correlated with
   log P — therefore rests on four database numbers and four binding rates. **k4b hold-out guard #4 refuses
   "any shipped matrix term that is a monotone function of log P", and nothing here changes that.** The
   correlation is also confounded four ways at once: functional group (pyrazine / furanthiol / alkenal /
   thiazole), molecular flexibility (the paper's own competing explanation is **steric hindrance**, not log P),
   covalent reactivity, and **dose** (Flags 4). **No single-variable conclusion can be drawn from four
   compounds.**
4. **The four compounds are not dosed at a common concentration.** The thiazole is at **0.384 mmol/L** and the
   other three at **0.066-0.106 mmol/L (mine)** — a **3.6-5.8x** difference. Binding rate is an
   occupancy-dependent quantity, so four rates measured at four occupancies are not strictly comparable, and
   the compound with the second-highest binding rate is also the one dosed 4-6x higher. (The direction happens
   to work against the paper's ordering rather than for it — saturation would *lower* the thiazole's fractional
   binding — so the qualitative ordering survives; but the numbers are not on a common basis.) The paper never
   explains why the thiazole was dosed 5.5x higher in mass terms.
5. **The Klotz numbers have no molar basis and an unreported Hill coefficient.** Eq. 7 divides by Cp in mol/L
   and **no pea protein molar mass appears anywhere in the paper** — identical to Bi 2022's Flags 1. And Eq. 6
   is written with an exponent h, "the Hill coefficient describing the cooperative character between a single
   binding site and the flavour ligand", with the note that h = 1 recovers the Klotz equation. **The fitted h
   is never printed.** If h ≠ 1 then K is not in L/mol and nK is not an affinity in the usual sense; if h = 1
   it should be said. Beyond that, **n and K individually are not significantly different between the control
   and the heat-treated group** on the paper's own letters, with standard errors of 64-77 % of the values —
   **only nK is usable, and only as a ratio.**
6. **The equilibrium is at 25 °C and the sampling is at 50 °C.** After the heat treatment and the ice-water
   quench, the vial is shaken 12 h at 25 °C "使之达到平衡状态" (to reach the equilibrium state) — and then, for
   the SPME step, **equilibrated at 50 °C for 10 min and extracted for 30 more minutes at that temperature**
   (§1.3.2). **Forty minutes at 50 °C is ample time for a reversible binding equilibrium to re-establish
   itself at 50 °C**, so it is not obvious which temperature the reported partition actually describes. The
   paper does not address this. It matters less for a within-run ratio (both the sample and the blank see the
   same 50 °C) than for an absolute value, but it does mean **`temperature_c` for a Xu row is genuinely
   uncertain between 25 and 50 °C, and 25 is the conservative choice**.
7. **The correlations are on five points and one of them is r = 1.00.** The Spearman coefficients quoted in
   §2.3.2 — r = 1.00 with heat-treatment temperature, 0.90 with time, −0.96 with particle size, −0.89 with
   surface hydrophobicity — are computed across the **five** conditions of a treatment series (control plus
   four levels). A Spearman r of exactly 1.00 on five points means only that the ranks are monotone, which for
   a four-level temperature ladder is a weak claim, and n = 5 gives a two-sided p of 0.017 at best. **These are
   descriptive, not inferential**, and the paper's mechanistic argument (binding sites exposed on the
   hydrophobic surface) leans on them. Note also the internal tension: the pyrazine binding rate correlates at
   r = 1.00 with temperature across 80-120 °C, but **surface hydrophobicity PEAKS at 100 °C and falls at
   120 °C** (Fig. 3A) — so binding and surface hydrophobicity **diverge at exactly the highest temperature**,
   which the paper acknowledges only for the thiazole.
8. **The paper's quotation of Guo 2019 is loose.** §2.4.2 attributes to Guo "high-affinity primary sites
   n = 1-2, K = 14 000-20 000 L/mol; low-affinity secondary sites n = 1-11, K = 100-1 500 L/mol". Guo 2019's
   Table 2 (see `guo2019_extraction.md`) actually gives primary **n = 0.8-1.4** and **K = 14 000-40 000** and
   secondary **n = 6.0-11.0** and **K = 290-1 500**. **The upper bound on the primary K is understated by 2x
   and the lower bound on the secondary n by 6x.** The comparison Xu draws — that its pea sites are "more
   numerous and weaker" than Guo's soy sites, hence surface-located — survives the correction, but do not
   propagate the quoted ranges; go to Guo's table.
9. **Table 2's ΔH is a van 't Hoff enthalpy, not an activation energy.** −86.68 and −60.83 kJ/mol come from
   fitting **ln Ka against 1/T** over 298/303/310 K (Eq. 4). They are the enthalpies of a binding equilibrium
   measured by fluorescence quenching. **They are NOT E_a**, they must not be put in `matrix_sites.py`'s
   `ea_band_kj_mol` (which carries 15-20 kJ/mol from real rate measurements), and the magnitudes make the
   error inviting: −86.68 kJ/mol would look like a very plausible activation energy for a Maillard step. **An
   activation energy cannot be negative; this is.** There is no rate and no time axis anywhere in this paper
   apart from the 5-20 min heat-treatment ladder, which is a treatment duration and not a kinetic measurement.
   Separately, the paper **prints the gas constant as "8.314 kJ/(mol·K)"** (p. 30), which is wrong by 1000x —
   the arithmetic used J/(mol·K) and this dossier verifies it (section 3 item 4), but a reader copying the
   printed unit would be badly misled.
10. **The van 't Hoff fit is three points over 12 K, and its R^2 is the weakest in the paper.** The paper prints
    **R^2 > 0.93** for the van 't Hoff fits against **> 0.98** for Stern-Volmer and **> 0.99** for the modified
    form. Twelve kelvin is a narrow window for an enthalpy, and the ΔH difference between the two groups
    (86.68 vs 60.83 kJ/mol, a 30 % change) carries no uncertainty at all — **neither ΔH nor ΔS is printed with
    an error bar.**
11. **The fluorescence and headspace experiments are at different loadings and disagree by 12-13x.** The
    quenching runs at **1 mg/mL** and the headspace at **5 mg/mL**, a 5x difference on top of the method
    difference; the resulting constants differ by **12.0x and 12.7x (mine, section 3 item 5)**. Bi 2022 — same
    laboratory, same protein, same buffer — showed 5-56x on the same comparison. **Do not present a Xu number
    without its method.** The one encouraging note is that both methods agree on the direction and roughly the
    size of the preheat effect (1.66x by Ka, 1.75x by nK).
12. **Almost every binding rate is figure-only.** Figs. 1 and 2 hold forty bars between them and **five numbers
    are printed**, all for two of the four compounds. So there is **no printed binding rate at 80 °C or 90 °C
    for any compound**, none at 15 or 20 min, none for the furanthiol or the decadienal beyond ">98 %", and no
    numeric value for the thiazole at 120 °C. The three reconstructions in section 3 (2-methylpyrazine at
    120 °C and at 15 min, the thiazole at 5 min) come from the text's "about X percentage points" phrases and
    **are estimates, not printed values**. Requesting the underlying table would roughly quadruple what this
    paper contributes.
13. **The protein is characterised by three bulk assays and nothing chemical.** Surface hydrophobicity,
    particle size and zeta potential — no SDS-PAGE, no DSC, no thiol, no disulfide, no free amine, no
    solubility, and no molar mass. `protein_matrices.yml`'s `pea_isolate` densities come from Gao 2020, Xiao
    2024, Chen 2022, Shen 2022 and Chihi 2016 on other preparations, so pairing a Xu binding ratio with those
    site densities is a cross-preparation pairing — **and it matters more here than usual, because the covalent
    claim for the furanthiol depends on the disulfide density and the claim for the decadienal on the amine
    density, neither of which was measured on this isolate.**
14. **Methanol is present at roughly 1 % v/v and is never accounted for (mine, section 3 item 8).** The stocks
    are 0.01 mol/L in methanol; the paper does not state the final co-solvent fraction, runs no methanol-only
    control, and never mentions it again — the same defect as Bi 2022 (~1.25 %) and Guo 2019 (propylene glycol).
15. **The isolate is at 5 g/L here and 10 g/L in Bi 2022, from the same laboratory.** So the two papers'
    per-gram constants are on nominally the same basis but were measured at a 2x different loading, and the
    linearity of K_g in loading is assumed, not demonstrated, in either. A `pea_protein_1pct` medium string
    already exists in `REVERSIBLE_BINDING` for the Bi rows; **a Xu row is at 0.5 % and needs its own string.**
16. **What this paper does NOT contain**: any binding measurement above 37 °C (fluorescence) or 50 °C (SPME
    sampling); any pea protein molar mass; any fitted Hill coefficient; any rate constant or activation energy;
    any adduct identification; any reversibility control; any binding constant for three of its four compounds;
    any pH other than 7.6; any ionic-strength series; any protein-free heated blank described as such; any
    sensory measurement or odour threshold; any supplementary material.
17. **What to request from the authors** (Wu Jihong, wjhcau@hotmail.com — the same corresponding author as
    Bi 2022, so one message could resolve both papers): (i) **whether the protein-free blank was heat-treated
    alongside the sample**, which decides how much of the temperature series is real binding; (ii) **a 2-methylpyrazine
    mass balance for the 120 °C arm**, since the compound is a Maillard product; (iii) **the pea protein molar
    mass used for Cp** — this would unlock the Klotz constants of both this paper and Bi 2022 at a stroke;
    (iv) the fitted Hill coefficient h; (v) the forty binding rates behind Figs. 1 and 2 as a table; (vi) the
    surface-hydrophobicity, particle-size and zeta-potential values as numbers; (vii) the final methanol
    fraction; (viii) an adduct search on the furanthiol and the decadienal, which would convert the paper's
    cited covalent claim into a measured one.
