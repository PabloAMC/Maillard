# Sun 2025 — EXTRACTION (pea protein isolate 18 g/L in 0.01 M PBS pH 7.2; dimethyl disulfide, dimethyl trisulfide and lenthionine by HS-SPME-GC-MS depletion at 50 C, fluorescence quenching at 298/304/310 K with a two-constant modified Stern-Volmer fit, van 't Hoff thermodynamics, bond-breaking reagents, zeta/particle size, UV, CD, AFM, E-nose and CDOCKER docking)

### THE FIRST SULFUR ROW THE BINDING TABLE COULD EVER HAVE: `REVERSIBLE_BINDING` in `src/kinetic_core/parameters_matrix.py` holds 21 constants over five media and **not one is a thiol, a disulfide, a trisulfide or any sulfur compound at all**, while the sulfur lane is the model's central output — and this paper measures pea protein isolate against **dimethyl disulfide, dimethyl trisulfide and lenthionine**, prints a headspace-depletion ladder, a fluorescence binding constant ladder, and a bond-breaking release test that says **DMDS is released by urea, i.e. bound non-covalently**. But the two ladders disagree by **1100x** on the same three compounds in the same paper (retention spans 1.57x, Ka spans 1755x, mine), and the paper prints **no partition coefficient and no water leg**, so nothing here converts to `K_g` without an assumption I have to declare.

**Source on disk:** `data/articles/Sun2025.pdf` (13 pp., Food Hydrocolloids 166 (2025) 111326).
Read from the `pdftotext -layout` text layer
(`scratchpad/Sun2025.txt`); **Tables 1 and 2 came through clean** and are re-typed in full below.
**Tables S1 and S2 are supplementary and are NOT on disk** — S1 holds the PEN3 sensor key, S2 holds
the CDOCKER binding affinities (whose three values ARE quoted verbatim in §3.8 and are typed here).
**Fig. 1A (the retention-vs-concentration curve), Fig. 1B–D (zeta, particle size, PDI), Fig. 2A–I
(UV, second-derivative, fluorescence, CD, secondary structure), Fig. 3 (AFM with the Rq values),
Fig. 4A–C (E-nose radar, correlation matrix, bond-breaking peak areas), Fig. 5 (docking) and
Fig. S1 (the modified Stern-Volmer and double-log plots) are images**: every retention percentage
other than the six span endpoints quoted in §3.1, every Rq, every zeta potential, every particle
size and every bond-breaking peak area is **figure-only**. Repo status before this dossier: Sun 2025
is cited nowhere in `src/kinetic_core/parameters_matrix.py`, `src/kinetic_core/matrix_sites.py`,
`src/kinetic_core/matrix_oav.py` or `data/species/protein_matrices.yml`, and has no extraction
dossier.

## 0. Identity

| field | value |
|---|---|
| Title | "Interactions of pea protein with three sulfur-containing flavor compounds: Insights into molecule structural, non-covalent, and binding mechanisms" |
| Authors | Hailan Sun, Jingyi Liang, Yirong Qian, Xiao Chen, Liyan Zhao (corresponding, zhlychen@njau.edu.cn) — College of Food Science and Technology, Nanjing Agricultural University, Nanjing 210095, China |
| Venue | Food Hydrocolloids **166** (2025), Article **111326**. Received 21 January 2025; revised 5 March 2025; accepted 5 March 2025; online 7 March 2025 |
| DOI | **`https://doi.org/10.1016/j.foodhyd.2025.111326`** — printed exactly so in the footer of p. 1 and again under "Appendix A. Supplementary data" on p. 12 |
| Funding | NSFC Youth Program 32402279; NSFC 32472477; PAPD (Priority Academic Program Development of Jiangsu Higher Education Institutions) |
| The three ligands | **DMD = dimethyl disulfide** (linear, analytically pure, J&K Scientific, PubChem ID 12,232); **DMT = dimethyl trisulfide** (linear, analytically pure, J&K, PubChem ID 19,310); **LEN = lenthionine** = 1,2,3,5,6-pentathiepane, printed in this paper as "1,2,3,5,6-pentathiolane" (cyclic, purity 98 %, Aladdin, PubChem ID 67,521) |
| Protein | **Commercial pea protein isolate, purity 90 %, Yuanye Biotechnology Co. Ltd (Shanghai)**. No further characterisation: no N x 6.25 assay of their own, no SDS-PAGE, no thiol, no free amine, no molar mass. Docking used pea 11S, **PDB 3KSC** — the same receptor Bi 2022 used |
| Naming | "Flavor retention rate (%)" = headspace **depletion** by Eq. 1 (see Flags 1 — the name is misleading); `KSV` = Stern-Volmer quenching constant; `KD`/`KS` = dynamic and static quenching constants from the modified Stern-Volmer Eq. 3; `Ka` = binding constant from the double-log Eq. 4 (Hill form), printed in units of `10^3/mol`; `n` = Hill coefficient (**not printed anywhere**); `H0` = ANS surface hydrophobicity index; `Rq` = AFM surface roughness |
| Companions on disk | `bi2022_extraction.md` (the same protein family, the same 3KSC receptor, the same bond-breaker method, the same group's method lineage — cited 6 times by this paper), `guo2020_extraction.md` (soy, cited as the 12 h/4 C equilibration source), `barallatperez2024_extraction.md` (lupin/commercial isolates, cited), `anantharamkrishnan2020b_extraction.md` (the DMDS adduct contradiction, cited by this paper as its opening animal-protein reference) |

## 1. Why it matters

**The sulfur gap in the matrix layer, stated precisely.** `REVERSIBLE_BINDING` in
`src/kinetic_core/parameters_matrix.py` currently carries 21 `MatrixParameter` rows across five
media — `skim_milk` (Meynier 2002), `caseinate_1pct` (Leksrisompong 2010), `soy_protein`
(Damodaran 1981 / Arai 1970), `beta_lactoglobulin` (Andriot 2000) and, since wave B26,
`pea_protein_1pct` (Bi 2022). Their compounds are esters, ketones, aldehydes, one lactone, one
alcohol and furaneol. **Not one is a thiol, a disulfide, a trisulfide, a thiophene or a thiazole.**
Meanwhile the repository's sulfur lane — 2-methyl-3-furanthiol, 2-furfurylthiol and their disulfides
— is the model's central output. Every OAV this model reports on its most important compound class
is therefore computed with a matrix shift of exactly 1.0, by absence of evidence rather than by
measurement. **This paper is the first corpus source that puts a plant protein against sulfur
volatiles and reports a number.** It does not close the gap (see the "but" below), and it is the
closest thing the corpus has.

**What it says about DMDS reversibility, which is the load-bearing question.**
`parameters_matrix.py` carries `COMPOUND_STRUCTURE["dimethyl_disulfide"]` with the binding class
`disulfide`, listed in the module's adduct-NEGATIVE half, and it carries the contradiction verbatim
in `SOURCE_CONTRADICTIONS["dimethyl_disulfide_adduct"]`: Anantharamkrishnan 2020b's Table 2 lists
DMDS as `no` (unreactive, one of the 32) while its own Results text on pp. 12-13 says DMDS **does**
form a covalent bond with beta-lactoglobulin and names the adduct **+46 Da, BLG-CysSSMe**. The
module resolves that conservatively — no covalent term for DMDS — and reports the contradiction on
every DMDS prediction. **Sun 2025 is directly relevant to that unresolved question and it comes down
on the non-covalent side, but only qualitatively and only on a different protein.** Its evidence:

1. The title, abstract and conclusion all say **non-covalent** ("the dominance of non-covalent
   forces", §4).
2. **§3.7, the bond-breaker experiment, is the only actual reversibility test, and DMDS passes it**:
   "Urea addition significantly increased DMD's peak area in the headspace (p < 0.05), confirming
   strong hydrophobic interaction with PPI." A compound that had left the headspace by forming a
   +46 Da covalent adduct with a cysteine would **not** come back when 4 M urea is added — urea
   unfolds, it does not cleave a disulfide interchange product. **That the DMDS peak area rises
   significantly on urea is a positive demonstration that a material fraction of the DMDS depletion
   was reversible.** The paper does not print the size of that fraction (Fig. 4C is an image), so
   the test is directional only.
3. The thermodynamics assign DMDS binding to **hydrogen bonds and van der Waals forces**
   (ΔH = −47.77 kJ/mol < 0 and ΔS = −114.49 J/mol·K < 0, Ross & Subramanian class 2), not to
   covalent chemistry.

**What it does NOT settle.** This is pea protein, not beta-lactoglobulin, and the contradiction in
`SOURCE_CONTRADICTIONS` is about beta-lactoglobulin. Sun performs **no mass spectrometry for
adducts**, searches for no +46 Da species, and reports no protein-side measurement at all — the
whole DMDS story is the disappearance of a GC peak and its partial return under urea. It therefore
**corroborates the conservative resolution already shipped (DMDS gets no covalent term) without
resolving the Anantharamkrishnan contradiction**, and it should be recorded in
`SOURCE_CONTRADICTIONS` as a corroborating-but-not-deciding third observation, on a different
protein, by a method that cannot see an adduct. Note the asymmetry the repository already encodes:
`ADDUCT_POSITIVE_CLASSES` contains `"trisulfide"`, so **DMTS already has a covalent gate open in this
model while DMDS does not** — and Sun treats DMD and DMT identically, as two linear thioethers,
with no covalent control on either. Sun's DMT numbers are therefore **more** contaminated by
irreversible chemistry than its DMD numbers, by the repository's own gate, and Sun's own
thermodynamics cannot tell the difference.

**The number the repository actually wants, and why this paper does not print it.** The stored form
is `K_g = (K_water/K_matrix − 1) / protein_g_per_L` in L/g protein, and every FIT row in the table
was built the same way: a **matrix leg and a water leg measured in the same run on the same
instrument**, so that the suspect absolute static-headspace scale cancels (Meynier is 6.24x low
against its own printed Henry's constants; Leksrisompong 6-17x low; both are shipped as ratios for
exactly that reason). **Sun runs a matched control** — "For the control group, the PBS solution was
used in place of the PPI solution, while all other conditions remained the same" (§2.2) — which is
structurally the right construction. But it reports the result only as a **depletion percentage**,
never as a partition coefficient, and the percentages themselves are in Fig. 1A. From the three
maxima quoted in the §3.1 text plus the printed statement that retention peaked at 0.5 mM, a `K_g`
can be reconstructed (§3 item 3, all mine) — **but it is SPME, not static headspace**, and SPME is a
competitive, fibre-capacity-limited, non-equilibrium sampling method whose depletion percentage is
not a partition coefficient (Flags 3). **The reconstructed constants are ~1.8-3.5e-2 L/g, i.e. an
order of magnitude BELOW Bi 2022's shipped pea hexanal row of 2.537e-1 L/g** — a sulfur-vs-aldehyde
contrast on the same protein species that is worth recording as a bound even if the absolute is not
shippable.

**Where it lands in the method taxonomy.** The registry's `method` field is first-class because k2
sec. B.3 measured a 35x aldehyde gap between headspace-depletion and dialysis determinations, and
`binding_constant_for` refuses to cross that boundary for aldehydes. Sun contributes **a third
method family the registry does not yet name: `hs_spme_depletion`.** It must not be filed as
`static_headspace_partition` (the Meynier/Leksrisompong/Bi-PRV method) and it must not be filed as
`headspace_depletion` (Andriot's equilibrium static-headspace depletion) without a note, because
SPME pre-concentrates onto a fibre with finite, compound-dependent, competitive capacity. And Sun's
fluorescence `Ka` is a **fourth object again** — a quenching-derived constant at 0.2 mg/mL, 90x
below the 18 g/L headspace loading. **This paper reproduces the Bi 2022 lesson on its own data and
more violently**: see §3 item 4, where the two ladders disagree by 1100x.

**What this paper does NOT give the repository**: any odour threshold measured by the authors (the
one threshold it prints, 0.008 ug/kg for DMDS, is **cited**); any partition coefficient; any water
leg as a number; any rate constant; any activation energy; any adduct mass; any pH other than 7.2;
any temperature above 50 C (and 50 C is only the SPME incubation, not the binding condition); any
2-methyl-3-furanthiol, 2-furfurylthiol or any other thiol; any heat-treated protein; any Hill
coefficient `n`.

## 2. Methods as they matter to a model

- **The pot.** A **20 mL headspace vial**, sealed, for the binding assay and the bond-breaker
  assay. 1.8 mL of PPI solution + 0.2 mL of SCFC stock = **2.0 mL liquid in 20 mL**, so
  beta ≈ 9 (mine, from the printed volumes). Fluorescence, UV, CD, zeta and AFM are in cuvettes /
  on mica, at far lower protein.
- **Protein loading — the number every per-gram constant divides by.** "The 2 % (w/v) solution of
  PPI was prepared by dissolving the protein in 0.01 M PBS (pH 7.2)"; 1.8 mL of that is diluted to
  2.0 mL, and the paper states the result explicitly: **"Each sample vial contained a final PPI
  concentration of 1.8 % (w/v)"** = **18 g/L of isolate**. The isolate is **90 % pure**, so
  **16.2 g/L of actual protein (mine)**. Both numbers are given in §4 below; the repository's
  convention on the existing rows is grams of the material named by the source, and `skim_milk`'s
  33.9 g/L is a protein figure, so **16.2 g/L is the like-for-like number and 18 g/L is the
  as-printed one**.
- **Buffer and pH.** **0.01 M PBS, pH 7.2**, throughout. Bi 2022's pea rows are at pH 7.6 in
  0.01 M potassium phosphate — close but not the same, and the registry stores pH as a first-class
  field, so these are two different rows and not a replication.
- **Methanol is present, and at a higher level than in Bi 2022.** "Stock solutions of three SCFCs
  (10 mM) were prepared in **methanol** and sonicated for 1 h"; 0.2 mL of that goes into 2.0 mL, so
  the binding assay carries **10 % v/v methanol (mine)** — eight times Bi 2022's ~1.25 %. For the
  zeta/particle-size work the stocks are 5 mM in methanol and 0.1 mL goes into 1.0 mL, again
  **10 % v/v (mine)**. **The E-nose control substitutes methanol for the flavour compounds**, which
  shows the authors knew the co-solvent mattered there; the binding assay's control is PBS-for-PPI
  and therefore carries the same 10 % methanol on both legs, so the methanol at least cancels
  between control and sample. It does not cancel against any other paper (Flags 5).
- **Concentration ladder.** Eight concentrations, **0.1 mM to 1.0 mM**, for the binding assay.
  0.1-0.5 mM for UV and fluorescence. **0.5 mM** for surface hydrophobicity, zeta/particle size, CD,
  AFM, E-nose and the bond-breakers — 0.5 mM is the paper's chosen "critical flavor concentration"
  because retention peaked there.
- **Temperature and time — and they are not one condition.** The binding assay: vortexed
  **40 min**, then **equilibrated at 4 C for 12 h** (after J. Guo et al. 2020), then **incubated at
  50 C for 15 min**, then SPME extraction **25 min at 50 C**. So the protein-flavour equilibrium is
  established cold (4 C, 12 h) and the headspace is sampled hot (50 C, 40 min total on the fibre).
  **There is no single temperature to put in the registry's `temperature_c` field** (Flags 2).
  Fluorescence is at **298 K / 304 K / 310 K**, 40 min, water bath. UV, CD, zeta, AFM at **25 C**.
  Bond-breakers: PPI + reagent stirred 1 h at 25 C, then + flavour, **equilibrated at 37 C for 2 h**.
- **Method 1 — HS-SPME-GC-MS depletion (the binding measurement).** Modified from Xi et al. 2023.
  **CAR/DVB/PDMS fibre, 50/30 um**, conditioned 250 C / 30 min. Incubate 50 C / 15 min, extract
  25 min, desorb **220 C / 3 min splitless**, 1 uL injection. **HP-5MS 30 m x 0.25 mm x 0.25 um**
  (note: a non-polar column, where Bi 2022 used HP-WAX). Helium 1.0 mL/min constant flow. Oven:
  40 C hold 3 min -> 130 C at 10 C/min hold 3 min -> 250 C at 6 C/min hold 5 min, 2 min post-run.
  MS: 70 eV, **scan 30-450 m/z** (full scan, not SIM), quadrupole 150 C, source 230 C. Instrument
  named as "a TSQ 1310 gas chromatograph ... equipped with a TSQ 9000 triple quadrupole" — the
  model numbers as printed (a Thermo TRACE 1310 GC with a TSQ 9000 MS is what that describes).
  **Eq. 1: Flavor retention rate (%) = (H0i − H1i)/H0i x 100 %**, H0i = peak area without protein,
  H1i = with protein. **This is a depletion fraction and nothing else.** It cannot separate
  reversible partitioning from irreversible capture; that separation is attempted only by the
  bond-breakers of §3.7, and only qualitatively.
- **Method 2 — fluorescence quenching (a different object).** Excitation fixed at **260 nm**
  (unusual — Bi 2022 and most protein work use 280-295 nm to select tryptophan; 260 nm excites
  more than Trp and sits on the paper's own UV absorption peak), emission **280-400 nm**, scan
  10 nm/s. **Inner-filter correction was done by subtracting SCFC-only control spectra** (after
  Guo 2024 and van de Weert & Stella 2011) — this is stated and is better practice than most of the
  corpus. Protein at **0.2 mg/mL** for the UV/fluorescence solutions prepared per §2.5.1, i.e.
  **90x below the 18 g/L headspace loading (mine)**. Three temperatures, 40 min equilibration.
  - Eq. 2, classical Stern-Volmer: `F0/F = 1 + Kq*tau0*[Q] = 1 + KSV[Q]`, with tau0 taken as 1e-8 s.
  - Eq. 3, **modified Stern-Volmer**: `F0/F = (1 + KD[Q])(1 + KS[Q])`, a second-order polynomial
    fitted to the upward curvature. **KD and KS are printed as identical to the last decimal in all
    nine rows of Table 2** (Flags 6 — this is the single most alarming thing in the paper).
  - Eq. 4, **double-log / Hill**: `lg((F0−F)/F) = lg Ka + n lg[Q]`, giving the binding constant `Ka`
    and Hill coefficient `n`. **`Ka` is tabulated; `n` is never printed.**
- **Method 3 — van 't Hoff thermodynamics.** Eq. 5 `ln Ka = −ΔH/(RT) + ΔS/R`, Eq. 6
  `ΔG = −RT ln Ka = ΔH − TΔS`, with ΔH and ΔS treated as constant over 298-310 K. **The printed ΔG
  column is ΔH − TΔS from the regression line, not −RT ln(printed Ka)** — I verified both (§3
  "Arithmetic", item 2); they agree for LEN and DMD and disagree by 0.7-1.5 kJ/mol for DMT.
- **Method 4 — surface hydrophobicity.** ANS, after Bi et al. 2022. Excitation 390 nm, emission
  470 nm, **slits 5 nm** (Bi used 1 nm). H0 = initial slope of fluorescence vs protein
  concentration (mg/mL). Reported as a percentage of the flavour-free control (= 100 %).
- **Method 5 — zeta potential and particle size.** Zetasizer Nano ZS90 (Malvern). PPI diluted to
  **0.2 mg/mL** in PBS; 0.9 mL + 0.1 mL of 5 mM SCFC stock; stirred 25 C / 40 min.
- **Method 6 — UV.** TU-1900 double-beam, 200-400 nm, 0.1 nm interval, 25 C; second derivative
  d2A/dlambda2 in Origin; the amplitude ratio **r = a/b** on the 280-300 nm peaks/troughs is the
  reported quantity. **No r value is printed as a number** — only its direction of change.
- **Method 7 — CD.** J-1500 (Jasco). 1.8 mL of **0.2 mg/mL** PPI + 0.2 mL of **5 mM** flavour (so
  0.5 mM final), 25 C / 40 min. **190-260 nm**, 1 nm interval, 50 nm/min, **methanol as the blank**.
  Secondary structure by the **Yang model**.
- **Method 8 — AFM.** Dimension Edge (Bruker), tapping mode at 25 C, **0.02 mg/mL protein**,
  0.5 mM flavour, 10 uL on freshly cleaved mica for 20 min, dried 30 C / 2 h, 3 um x 3 um scans,
  NanoScope Analysis. Rq values are printed **in the Fig. 3 panels only**.
- **Method 9 — E-nose.** PEN3 (Airsense), ten metal-oxide sensors. Final 1.8 % w/v PPI, 0.5 mM
  SCFC; incubated 50 C / 15 min then equilibrated 25 C / 30 min. Auto-clean 180 s, zero 5 s, prep
  5 s, measure 90 s, interval 10 s, flow 400 mL/min. Sensor key in Table S1 (not on disk).
  **An E-nose is not a sensory panel and its response values are not odour activity values.**
- **Method 10 — bond breakers (the reversibility test).** After Bi et al. 2022. PPI diluted to
  **2 mg/mL** in each of **80 % propylene glycol (PG)**, **4 mol/L urea** and **0.6 mol/L NaCl**;
  stirred 25 C / 1 h. Then 1.8 mL of treated PPI + 0.2 mL of 0.5 mM SCFC stock in a 20 mL vial,
  sealed, vortexed, **equilibrated 37 C for 2 h**, read by HS-SPME-GC-MS. Control = PPI in PBS with
  no reagent. Reading rule as printed: no significant change = no corresponding bond; a significant
  **increase** = a weakened interaction force releasing the compound. **Note the protein here is
  2 mg/mL, 9x below the 18 g/L of the binding assay** — the bond-breaker run is not at the binding
  assay's loading.
- **Method 11 — docking.** Discovery Studio 2016 (BIOVIA), **CDOCKER** protocol with the **CHARMm**
  force field, receptor **PDB 3KSC** (pea 11S) after removing waters, adding hydrogens and
  optimising side chains; ligands from PubChem. **This is force-field docking, NOT DFT** — it is
  nonetheless in silico and no number from it is a measurement.
- **Replication and statistics.** Triplicate; one-way ANOVA, p < 0.05; SPSS 20.0; Origin 2021.

## 3. Tables re-typed

### Table 1 (p. 5). "Decrease rate of PPI surface hydrophobicity index in the presence of DMD, DMT, and LEN, as well as the structural formula and hydrophobic constant (Log P) of three SCFCs."

| SCFCs | PPI (0 mM) | 0.5 mM | Structural formulas | Log P |
|---|---|---|---|---|
| (PPI, no flavour) | 100 % `[M]` | – | – | – |
| DMD | – | **99.36 %** `[M]` | *(structure drawing, image only)* | **1.77** `[C]` |
| DMT | – | **91.70 %** `[M]` | *(image only)* | **2.93** `[C]` |
| LEN | – | **87.99 %** `[M]` | *(image only)* | **4.23** `[C]` |

The title says "Decrease rate", but the tabulated quantity is the **remaining** H0 as a percentage
of the flavour-free control, which the text confirms: "the H0 of PPI decreased to 99.36 %, 91.70 %,
and 87.99 %". The **decrease** is therefore 0.64 / 8.30 / 12.01 percentage points (mine). The Log P
column is `[C]`: no method for measuring it is given anywhere in §2, and these are literature
octanol/water constants. The structural formulas are drawings and carry no number.

### Table 2 (p. 7). "Binding and thermodynamic parameters of PPI interacting with DMD, DMT, and LEN at different temperatures."

| SCFCs | T (K) | KD (L/mol) | KS (L/mol) | Ra | Ka (10^3/mol) | Rb | ΔH (kJ/mol) | ΔS (J/mol·K) | ΔG (kJ/mol) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DMD | 298 | 319.71 `[F]` | 319.71 `[F]` | 0.9918 `[F]` | 0.27 `[F]` | 0.9967 `[F]` | −47.77 `[F]` | −114.49 `[F]` | −13.65 `[F]` |
|  | 304 | 230.52 `[F]` | 230.52 `[F]` | 0.9993 `[F]` | 0.14 `[F]` | 0.9986 `[F]` | | | −12.96 `[F]` |
|  | 310 | 218.97 `[F]` | 218.97 `[F]` | 0.9998 `[F]` | 0.13 `[F]` | 0.9973 `[F]` | | | −12.27 `[F]` |
| DMT | 298 | 744.70 `[F]` | 744.70 `[F]` | 0.9976 `[F]` | 15.57 `[F]` | 0.9959 `[F]` | −131.74 `[F]` | −364.15 `[F]` | −23.22 `[F]` |
|  | 304 | 679.67 `[F]` | 679.67 `[F]` | 0.9946 `[F]` | 2.30 `[F]` | 0.9979 `[F]` | | | −21.03 `[F]` |
|  | 310 | 588.38 `[F]` | 588.38 `[F]` | 0.9941 `[F]` | 2.01 `[F]` | 0.9932 `[F]` | | | −18.85 `[F]` |
| LEN | 298 | 3558.45 `[F]` | 3558.45 `[F]` | 0.9989 `[F]` | 473.81 `[F]` | 0.9963 `[F]` | −74.66 `[F]` | −142.10 `[F]` | −32.32 `[F]` |
|  | 304 | 3452.78 `[F]` | 3452.78 `[F]` | 0.9900 `[F]` | 243.11 `[F]` | 0.9951 `[F]` | | | −31.47 `[F]` |
|  | 310 | 3259.44 `[F]` | 3259.44 `[F]` | 0.9969 `[F]` | 147.71 `[F]` | 0.9983 `[F]` | | | −30.61 `[F]` |

Printed note under the table: "Ra is the linear correlation coefficient of KD and KS, and Rb is the
linear correlation coefficient of Ka."

Every cell in this table is `[F]` — a fitted parameter or a fit statistic, not a directly measured
quantity. **`Ka` is in units of 10^3 L/mol**, so DMD at 298 K is **270 L/mol**, DMT **15 570 L/mol**,
LEN **473 810 L/mol**. **KD = KS exactly in all nine rows** (Flags 6). The ΔH/ΔS blocks are printed
once per compound and apply down all three temperatures by the constancy assumption of Eq. 5.

### Table S2 (supplementary, values quoted verbatim in §3.8, p. 8)

| SCFC | CDOCKER binding affinity, kcal/mol |
|---|---:|
| LEN | **−7.77** *(in silico)* |
| DMT | **−7.32** *(in silico)* |
| DMD | **−7.23** *(in silico)* |

Printed as: "the binding affinity of PPI to SCFCs was in the order of LEN (−7.77 kcal/mol) > DMT
(−7.32 kcal/mol) > DMD (−7.23 kcal/mol)". **Force-field docking (CDOCKER/CHARMm), not DFT, and not a
measurement.** Recorded here only because the paper leans on it; nothing in §4 ships it.

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| retention span, **LEN**, over 0.1-1.0 mM | **11.01 % to 38.89 %** `[M]` | §3.1, p. 5 |
| retention span, **DMT**, over 0.1-1.0 mM | **7.18 % to 27.13 %** `[M]` | §3.1, p. 5 |
| retention span, **DMD**, over 0.1-1.0 mM | **3.04 % to 24.73 %** `[M]` | §3.1, p. 5 |
| where the retention maximum sits | **"the retention capacity of PPI for all SCFCs peaked at 0.5 mM"**; rising 0.1->0.5 mM, falling above 0.5 mM | §3.1, p. 4 |
| the concentration chosen for everything downstream | **0.5 mM** ("the critical flavor concentration") | §3.1, p. 5 |
| surface hydrophobicity H0 at 0.5 mM, DMD / DMT / LEN | **99.36 % / 91.70 % / 87.99 %** of the flavour-free control `[M]` | §3.2, p. 5 (= Table 1) |
| alpha-helix content, blank PPI | **11.4 %** `[F]` (CD + Yang model) | §3.4.4, p. 7 |
| alpha-helix content with 0.5 mM DMD / DMT / LEN | **8.9 % / 8.8 % / 8.7 %** `[F]` | §3.4.4, p. 7 |
| beta-sheet | "dramatically raised (p < 0.05)" — **no value printed** | §3.4.4 |
| the Bi 2022 comparison the authors themselves draw | Bi reported pea alpha-helix **17.9 % -> 15.4 %** on 1 mM hexanal | §3.4.4, p. 7 `[C]` |
| **odour threshold of dimethyl disulfide** | **0.008 ug/kg** `[C]` — **CITED from J. Sun et al. 2021, not measured here** | §1 (Introduction), p. 1 |
| zeta potential | negative before and after; **significantly decreased in absolute value** by all three SCFCs (p < 0.05); **no significant difference BETWEEN** DMD, DMT and LEN (p > 0.05) — **no value printed** | §3.3, p. 5, Fig. 1B |
| particle size and PDI | both increased; **LEN > DMT > DMD** on mean particle size — **no value printed** | §3.3, Fig. 1C-D |
| AFM roughness Rq | increased, "especially following LEN addition" — **values are inside the Fig. 3 panels only** | §3.5, Fig. 3 |
| E-nose, blank PPI | **W5S** (nitrogen oxides) highest, then **W2S** (alcohols, aromatics, short-chain alkanes) and **W1S** (short-chain alkanes) — "contribute significantly to the beany flavor" | §3.6, Fig. 4A |
| E-nose, + SCFCs | **W1W** (sulfides, terpenes) and **W2W** (aromatics, organic sulfides) "markedly increased ... dominating the flavor profile"; **W2S decreased** after DMD and LEN | §3.6, Fig. 4A |
| correlation analysis | retention correlates **positively** (p < 0.05) with zeta potential, particle size, PDI, beta-sheet and Rq; **negatively** with alpha-helix and beta-turn; **positively** with W1W and W2W; **negatively** with W5C, W3C and W2S — **no coefficient printed** | §3.6, Fig. 4B |
| **bond breaker, DMD** | **urea SIGNIFICANTLY INCREASED DMD's peak area (p < 0.05)** -> hydrophobic interaction; PG slightly decreased it, **not significantly** -> hydrogen bonding also present | §3.7, p. 8, Fig. 4C |
| **bond breaker, DMT** | **PG SIGNIFICANTLY DECREASED DMT's peak area (p < 0.05)** -> strong hydrogen bonding; **urea SIGNIFICANTLY INCREASED it (p < 0.05)** -> hydrophobic interaction. "DMT-PPI interaction is primarily hydrogen bonding, followed by hydrophobic interactions" | §3.7, p. 8 |
| **bond breaker, LEN** | **urea SIGNIFICANTLY INCREASED LEN's peak area (p < 0.05)** -> strong hydrophobic interaction; PG slightly increased it, **not significant** | §3.7, p. 8 |
| **bond breaker, NaCl, all three** | **NO significant change (p > 0.05)** -> minimal electrostatic interaction | §3.7, p. 8 |
| **every bond-breaker peak area** | **figure-only (Fig. 4C)** — not one is printed as a number | §3.7 |
| docking, H-bond to **SER233** backbone hydrogen and the S atoms of DMD / DMT / LEN | **2.43 Å / 2.85 Å / 2.23 Å** *(in silico)* | §3.8, p. 8-10 |
| docking, LEN pi-sulfur bond to **PHE90** benzene ring | **5.86 Å** *(in silico)*, named "the critical factor for the most vital adsorption capacity of PPI to LEN" | §3.8 |
| docking, **ARG469** van der Waals with LEN; **ILE235** (printed "IIE235") alkyl interaction with both straight-chain thioethers | no distance printed | §3.8 |
| UV | absorption peak of PPI near **260 nm** rises with SCFC concentration; **LEN causes a slight red shift and induces a NEW peak near 328 nm** ("may be caused by the combination of protein and LEN to form a new complex, and the specific mechanism needs to be further studied") | §3.4.1, p. 6 |
| second-derivative amplitude ratio r = a/b | **decreased** with DMT and LEN — no value printed | §3.4.1 |
| fluorescence | PPI peak emission **300-320 nm**; quenching rises with SCFC concentration; **LEN quenches most strongly** | §3.4.2, p. 6 |
| quenching mechanism | modified Stern-Volmer plots curve **upward**; KD and KS both **decrease** with temperature -> **static quenching dominant**; curvature strongest for LEN | §3.4.2, p. 6, Fig. S1 |
| interaction assignment from thermodynamics | **ΔH < 0 and ΔS < 0 for all three** -> hydrogen bonds and van der Waals forces (Ross & Subramanian class 2) | §3.4.3, p. 7 |

### Arithmetic on the printed constants (all mine)

**1. The two protein loadings.** The vial holds **1.8 % w/v isolate = 18 g/L**, and the isolate is
**90 % pure**, so **16.2 g/L of protein (mine)**. Liquid 2.0 mL in a 20 mL vial gives
**beta = 9 (mine)**. Methanol: 0.2 mL of a methanolic stock in 2.0 mL = **10 % v/v (mine)**.

**2. Table 2's internal consistency.** Two independent checks:

*(a) ΔG against ΔH − TΔS, using the printed ΔH and ΔS:*

| compound | 298 K (mine) | 304 K (mine) | 310 K (mine) | printed |
|---|---:|---:|---:|---|
| DMD | −13.65 | −12.97 | −12.28 | −13.65 / −12.96 / −12.27 ✓ |
| DMT | −23.22 | −21.04 | −18.85 | −23.22 / −21.03 / −18.85 ✓ |
| LEN | −32.31 | −31.46 | −30.61 | −32.32 / −31.47 / −30.61 ✓ |

**The ΔG column is ΔH − TΔS to within rounding in all nine rows.** So ΔG was computed from the
van 't Hoff line, not from the individual Ka values.

*(b) ΔG against −RT ln(printed Ka):*

| compound | 298 K (mine) | 304 K (mine) | 310 K (mine) | printed ΔG | worst gap |
|---|---:|---:|---:|---|---:|
| DMD | −13.87 | −12.49 | −12.55 | −13.65 / −12.96 / −12.27 | **0.47 kJ/mol** |
| DMT | −23.92 | −19.56 | −19.60 | −23.22 / −21.03 / −18.85 | **1.47 kJ/mol** |
| LEN | −32.38 | −31.34 | −30.68 | −32.32 / −31.47 / −30.61 | 0.13 kJ/mol |

**LEN closes; DMD is loose; DMT does not close.** Eq. 6 asserts `ΔG = −RT ln Ka = ΔH − TΔS`, and for
DMT the two halves of that identity differ by up to 1.47 kJ/mol, which is a **1.8x discrepancy in
the implied Ka (mine)**. The cause is visible in the Ka column: DMT's Ka falls **15 570 -> 2 300**
over six kelvin (**6.8x, mine**) and then only **2 300 -> 2 010** over the next six (**1.14x,
mine**). That is not a van 't Hoff line, it is one outlier and two clustered points, and the ΔH of
−131.74 kJ/mol is dominated by the 298 K value. **DMT's thermodynamics are the least trustworthy
block in the table** (Flags 7).

*(c) recovering ΔH from the Ka endpoints (mine, to check the regression):* slope
= Δ(ln Ka)/Δ(1/T) over 298-310 K, times −R. DMD: **−46.7** against a printed −47.77. DMT: **−131.0**
against −131.74. LEN: **−74.6** against −74.66. **All three reproduce to within 2 %**, confirming
that ΔH is a two-to-three-point van 't Hoff fit on the Ka column and inherits every problem the Ka
column has.

**3. Reconstructing the registry's `K_g` form (mine — and see Flags 3 and 4 before using it).**
The registry stores `K_g = (K_water/K_matrix − 1) / protein_g_per_L`, where K is an AIR/matrix
partition coefficient. Sun's Eq. 1 gives the depletion fraction `R = (H0 − H1)/H0` at matched
volumes with a matched PBS control, so the headspace peak-area ratio is `H1/H0 = 1 − R` and, **under
the assumption that SPME fibre uptake is proportional to gas-phase concentration and that both legs
sample the same phase ratio**, `K_water/K_matrix = 1/(1 − R)`. Taking the **0.5 mM** point, which is
the maximum of each ladder and therefore equals the upper span endpoint printed in §3.1:

| compound | R at 0.5 mM | K_water/K_matrix (mine) | K_g at 18 g/L isolate, L/g (mine) | K_g at 16.2 g/L protein, L/g (mine) |
|---|---:|---:|---:|---:|
| DMD | 24.73 % | 1.3286 | **1.83e-2** | **2.03e-2** |
| DMT | 27.13 % | 1.3723 | **2.07e-2** | **2.30e-2** |
| LEN | 38.89 % | 1.6364 | **3.54e-2** | **3.93e-2** |

**These three numbers are the whole reason to read this paper, and they are not shippable as
printed** — they are SPME depletion, not a partition coefficient, and they rest on an equality
between fibre uptake ratio and partition ratio that the paper never establishes (Flags 3). What they
are good for is a **magnitude bound**: pea protein at 16-18 g/L removes a quarter to two-fifths of
these three sulfur volatiles from the headspace, which puts a per-gram constant in the
**1e-2 to 4e-2 L/g** band. For comparison, the shipped `kg_hexanal_pea` is **2.537e-1 L/g** — so
**pea protein's apparent grip on these sulfur compounds is 6-14x weaker than its grip on hexanal
(mine)**, across two papers, two methods and two pH values.

**4. The two ladders in this paper disagree by 1100x, which is worse than Bi 2022's 56x.** Same
laboratory, same commercial protein lot, same buffer, same three ligands:

| ladder | DMD | DMT | LEN | LEN/DMD span |
|---|---:|---:|---:|---:|
| headspace retention at 0.5 mM (%) | 24.73 | 27.13 | 38.89 | **1.57x (mine)** |
| reconstructed K_g (L/g, mine) | 1.83e-2 | 2.07e-2 | 3.54e-2 | **1.93x (mine)** |
| fluorescence KSV = KD = KS at 298 K (L/mol) | 319.71 | 744.70 | 3558.45 | **11.1x (mine)** |
| fluorescence **Ka** at 298 K (L/mol) | 270 | 15 570 | 473 810 | **1755x (mine)** |
| fluorescence **Ka** at 310 K (L/mol) | 130 | 2 010 | 147 710 | **1136x (mine)** |
| CDOCKER affinity (kcal/mol, in silico) | −7.23 | −7.32 | −7.77 | 0.54 kcal/mol |

**The order is the same in every row — LEN > DMT > DMD — and the MAGNITUDE spans 1.57x on the
headspace and 1755x on the fluorescence, a factor of 1118 between the two determinations of the
same ladder (mine).** The paper repeatedly says the methods "align", and on rank order they do; on
magnitude they are not the same measurement at all. **A fluorescence-quenching-derived Ka is NOT a
headspace partition constant**, and this paper is now the corpus's second and most extreme
demonstration of it, after Bi 2022's 56x on hexanal. It should be cited alongside Bi 2022 in
`parameters_matrix.py`'s §2 header comment and in the notes on `binding_constant_for`.

**5. Within-compound temperature contrast (mine).** Ka falls with temperature for all three:
DMD 270 -> 130 (**2.08x over 12 K**), DMT 15 570 -> 2 010 (**7.75x**), LEN 473 810 -> 147 710
(**3.21x**). KSV falls too: DMD **1.46x**, DMT **1.27x**, LEN **1.09x**. **Ka and KSV disagree on
how temperature-sensitive the same interaction is, by up to 6x (DMT, mine)**, from the same
fluorescence spectra fitted two different ways.

**6. Structure response (mine).** Alpha-helix falls **11.4 % -> 8.9 / 8.8 / 8.7 %**, i.e. absolute
losses of **2.5 / 2.6 / 2.7 percentage points** and relative losses of **21.9 % / 22.8 % / 23.7 %**.
The three are nearly identical — **the conformational response does NOT resolve the ladder**, even
though H0 and particle size do. Surface hydrophobicity falls **0.64 / 8.30 / 12.01 percentage
points**, a **19x span from DMD to LEN (mine)**, which is the only structural readout that tracks
the retention ladder with any dynamic range.

**7. Log P against the ladder — and why it must not be shipped.** Log P 1.77 / 2.93 / 4.23 against
retention 24.73 / 27.13 / 38.89 % is monotone, and the paper builds its whole structural argument on
it. **`parameters_matrix.py` refuses "any shipped matrix term that is a monotone function of log P"
(k4b hold-out guard #4).** Three points, all monotone in three other variables too (sulfur count,
ring vs chain, molecular volume), license nothing. Recorded so a later reader does not rediscover it
and mistake it for evidence.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** `dimethyl_disulfide` is keyed (line 675) and has a
`COMPOUND_STRUCTURE` entry in `parameters_matrix.py` (binding class `disulfide`, carrying the
Anantharamkrishnan contradiction in its note). `dimethyl_trisulfide` is keyed (line 688) and has
**no** `COMPOUND_STRUCTURE` entry, though the class name `trisulfide` already appears in
`ADDUCT_POSITIVE_CLASSES`. **`lenthionine` is NOT in `data/keys/compounds.yml` at all** and would be
the first cyclic polysulfide anywhere in the repository.

Every row below shares: **commercial pea protein isolate (Yuanye, purity 90 %), 0.01 M PBS pH 7.2,
10 % v/v methanol, triplicate.** The headspace rows are at **18 g/L isolate = 16.2 g/L protein**,
vortexed 40 min, **equilibrated 4 C / 12 h**, then **sampled by SPME at 50 C**. The fluorescence rows
are at **0.2 mg/mL**, 40 min, at the stated temperature.

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| headspace depletion ("flavor retention rate"), **DMD** | **24.73** (max; span 3.04-24.73 over 0.1-1.0 mM) | % of headspace peak area removed | 0.5 mM, 18 g/L isolate, pH 7.2, 4 C 12 h then SPME at 50 C | §3.1, p. 5 (span endpoints; the eight-point curve is Fig. 1A) | **level_only** — a depletion fraction, not a constant |
| headspace depletion, **DMT** | **27.13** (span 7.18-27.13) | % | as above | §3.1, p. 5 | **level_only** |
| headspace depletion, **LEN** | **38.89** (span 11.01-38.89) | % | as above | §3.1, p. 5 | **level_only** |
| the maximum sits at 0.5 mM | "peaked at 0.5 mM"; rising below, falling above | — | as above | §3.1, p. 4 | **measured_bound** — this is what licenses reading the upper span endpoint as the 0.5 mM value |
| **reconstructed per-gram binding constant, DMD** | **2.03e-2** (1.83e-2 on the as-printed 18 g/L isolate basis) | L/g protein | 0.5 mM, 16.2 g/L protein, pH 7.2, SPME depletion | `(1/(1−0.2473) − 1)/16.2` from §3.1 **(mine)** | **derived_assumption** — SPME uptake assumed proportional to gas concentration (Flags 3) |
| **reconstructed per-gram binding constant, DMT** | **2.30e-2** (1.83e-2 -> 2.07e-2 on the 18 g/L basis) | L/g protein | as above | **(mine)** | **derived_assumption** |
| **reconstructed per-gram binding constant, LEN** | **3.93e-2** (3.54e-2 on the 18 g/L basis) | L/g protein | as above | **(mine)** | **derived_assumption** |
| **magnitude bound: any pea x sulfur per-gram constant from this paper** | **1e-2 to 4e-2** | L/g protein | 16-18 g/L, pH 7.2, SPME | §3 item 3 **(mine)** | **measured_bound** — the defensible form of the three rows above |
| **sulfur-vs-aldehyde contrast on pea protein** | shipped `kg_hexanal_pea` 2.537e-1 L/g is **6.2-13.9x** the three reconstructed sulfur values | — | cross-study (Bi 2022 pH 7.6 PRV vs Sun 2025 pH 7.2 SPME) | **(mine)** | **derived_assumption** — cross-study AND cross-method; a magnitude statement only |
| binding constant Ka, **DMD** | **270 / 140 / 130** (0.27 / 0.14 / 0.13 x10^3) | L/mol | 0.2 mg/mL, 298 / 304 / 310 K, **fluorescence quenching**, double-log Hill fit | Table 2, p. 7 | **binding_constant** — **fluorescence-derived; must NOT be pooled with any headspace value** |
| binding constant Ka, **DMT** | **15 570 / 2 300 / 2 010** | L/mol | as above | Table 2, p. 7 | **binding_constant**, fluorescence-derived. Note the 6.8x drop over the first 6 K (Flags 7) |
| binding constant Ka, **LEN** | **473 810 / 243 110 / 147 710** | L/mol | as above | Table 2, p. 7 | **binding_constant**, fluorescence-derived |
| Hill coefficient n (the companion of Ka in Eq. 4) | — | — | — | never printed | **absent** — Ka in L/mol is only interpretable with n, and n is missing (Flags 8) |
| Stern-Volmer / static / dynamic quenching constants KD = KS, DMD | **319.71 / 230.52 / 218.97** | L/mol | 0.2 mg/mL, 298 / 304 / 310 K, fluorescence | Table 2, p. 7 | **binding_constant** (quenching family) — **KD and KS are printed identical, Flags 6** |
| KD = KS, DMT | **744.70 / 679.67 / 588.38** | L/mol | as above | Table 2, p. 7 | as above |
| KD = KS, LEN | **3558.45 / 3452.78 / 3259.44** | L/mol | as above | Table 2, p. 7 | as above |
| ΔH / ΔS / ΔG(310 K), **DMD** | **−47.77 / −114.49 / −12.27** | kJ/mol, J/mol·K, kJ/mol | van 't Hoff on Ka over 298-310 K, 0.2 mg/mL | Table 2, p. 7 | **binding_constant** (thermodynamic companion, fluorescence-derived). **ΔH is a van 't Hoff enthalpy, NOT an activation energy** |
| ΔH / ΔS / ΔG(310 K), **DMT** | **−131.74 / −364.15 / −18.85** | as above | as above | Table 2, p. 7 | as above, and the least trustworthy block (Flags 7) |
| ΔH / ΔS / ΔG(310 K), **LEN** | **−74.66 / −142.10 / −30.61** | as above | as above | Table 2, p. 7 | as above |
| **reversibility of DMDS binding to pea protein** | **4 M urea significantly INCREASES the DMD headspace peak area (p < 0.05)** | qualitative | 2 mg/mL PPI, 0.5 mM DMD, 37 C 2 h | §3.7, p. 8 (Fig. 4C carries the areas) | **measured_bound** — a positive, directional demonstration that a material share of DMDS depletion is reversible. **The share is NOT quantified** |
| reversibility, DMT | **PG significantly DECREASES** DMT's peak area (H-bonding); **urea significantly INCREASES** it (hydrophobic) | qualitative | as above | §3.7, p. 8 | **measured_bound** |
| reversibility, LEN | **urea significantly INCREASES** LEN's peak area; PG's increase is not significant | qualitative | as above | §3.7, p. 8 | **measured_bound** |
| electrostatic contribution, all three | **0.6 M NaCl: no significant change (p > 0.05)** | qualitative | as above | §3.7, p. 8 | **structural_gate** — electrostatics excluded as a binding force for all three sulfur compounds on pea protein. Corroborated independently by the zeta result (§3.3: no significant difference between the three, p > 0.05) |
| surface hydrophobicity H0 remaining at 0.5 mM, DMD / DMT / LEN | **99.36 / 91.70 / 87.99** | % of flavour-free control | ANS, ex 390 / em 470 nm, slits 5 nm | Table 1, p. 5 | **level_only** |
| alpha-helix, blank / DMD / DMT / LEN | **11.4 / 8.9 / 8.8 / 8.7** | % | 0.2 mg/mL, 0.5 mM, 25 C, CD 190-260 nm, Yang model | §3.4.4, p. 7 | **level_only** — and it does NOT resolve the ladder (§3 item 6, mine) |
| Log P, DMD / DMT / LEN | **1.77 / 2.93 / 4.23** | — | literature octanol/water | Table 1, p. 5 | **level_only**, and `[C]` — **NOT to be used**: k4b hold-out guard #4 forbids any shipped log-P-monotone term |
| odour threshold, dimethyl disulfide | **0.008** | ug/kg | medium not stated by this paper | §1, p. 1 — **CITED from J. Sun et al. 2021** | **threshold `[C]`** — a second-hand citation with no medium, no panel and no criterion. **Do not enter it as a threshold; chase the primary source if the value is wanted** |
| within-study ladder ratio, headspace retention, LEN/DMD | **1.57** | — | 0.5 mM, 16.2 g/L, pH 7.2, SPME | 38.89/24.73 **(mine)** | **within_study_ratio** |
| within-study ladder ratio, fluorescence Ka at 298 K, LEN/DMD | **1755** | — | 0.2 mg/mL, 298 K | 473 810/270 **(mine)** | **within_study_ratio** |
| **method disagreement, headspace ladder vs fluorescence Ka ladder** | **1118x** on the same three compounds in the same paper | — | same lab, same lot, same buffer, same ligands | 1755/1.57 **(mine)** | **within_study_ratio** — the corpus's most extreme single-paper demonstration of the `method` boundary, ahead of Bi 2022's 56x |
| ladder ratio, DMT/DMD, headspace | **1.097** | — | 0.5 mM | 27.13/24.73 **(mine)** | **within_study_ratio** — the disulfide/trisulfide step is only 10 % on the headspace and **57.7x** on Ka **(mine)** |
| CDOCKER binding affinity, LEN / DMT / DMD | **−7.77 / −7.32 / −7.23** | kcal/mol | in silico, CHARMm force field, PDB 3KSC | Table S2 via §3.8, p. 8 | **derived_assumption** — in silico, **not DFT**, not a measurement, not to be shipped |
| docking H-bond distances, SER233 to DMD / DMT / LEN | **2.43 / 2.85 / 2.23** | Å | in silico | §3.8, p. 10 | **derived_assumption** |
| E-nose sensor responses, zeta potentials, particle sizes, PDI, Rq, all bond-breaker peak areas, all correlation coefficients, the second-derivative r values, the eight-point retention curve | — | — | — | Figs. 1A-D, 2A-I, 3, 4A-C, S1 | **figure_only** |

### Can these be put on the same basis as the shipped binding constants?

**(a) The reconstructed K_g CAN be put in the registry's form arithmetically, and SHOULD NOT be
shipped as a `measured_ratio`.** Every FIT row in `REVERSIBLE_BINDING` is a within-run ratio of two
legs on one instrument, and Sun's matched PBS control is that construction. But all five existing
FIT sources measure a **partition coefficient** (static headspace, PRV, dialysis, gel filtration);
Sun measures **SPME fibre uptake**, which is competitive, capacity-limited and not proportional to
gas concentration in general. The honest filing is a **new `measured_bound` row per compound** with
`method = "hs_spme_depletion"` and the band 1e-2 to 4e-2 L/g, or nothing at all. **Do not create a
`MatrixParameter` with `method = "static_headspace_partition"` from this paper.**

**(b) There is no single temperature to file.** The registry's `temperature_c` is a scalar. Sun's
binding equilibrium is established at **4 C over 12 h** and read at **50 C over 40 min**; the
bond-breaker run is at **37 C / 2 h**; the fluorescence at **25 / 31 / 37 C**. Bi 2022's shipped pea
rows are a clean 37 C. **Filing Sun at any one temperature is a fabrication.** If a row is created,
`temperature_c` should be left `None` and the schedule written into `notes`.

**(c) pH 7.2 is not pH 7.6.** The shipped `pea_protein_1pct` loading carries pH 7.6 from Bi 2022, and
`MATRIX_LOADING` stores pH as a field. Sun needs a **second pea loading entry** —
`pea_protein_1pct8` or similar at 16.2-18 g/L, pH 7.2 — not a reuse of Bi's.

**(d) The Ka values cannot be put on a per-gram basis at all.** They are in L/mol and come from a
Hill fit (Eq. 4) whose **exponent n is never printed**, so they are not even dimensionally
interpretable as an ordinary association constant unless n = 1, which for a 473 810 L/mol value on a
473 810-fold ladder is not credible. There is no molar mass for pea protein in this paper either.
**No per-gram number can be derived from Table 2. Ship none of it.**

**(e) Nothing here goes to `matrix_sites.py`.** That module wants second-order rate constants in
M^-1 s^-1 and activation energies in kJ/mol. This paper has neither. Its ΔH values are van 't Hoff
enthalpies of an equilibrium; reading −131.74 kJ/mol as an Ea would be a category error, and the
sign would make the channel run backwards.

**(f) What SHOULD change in the repository today**, on this evidence:
1. `SOURCE_CONTRADICTIONS["dimethyl_disulfide_adduct"]` gains a third observation: Sun 2025 §3.7
   shows urea-releasable DMDS binding on **pea** protein, corroborating the conservative
   no-covalent-term resolution **without** resolving the beta-lactoglobulin contradiction, and by a
   method (GC peak area) that cannot see a +46 Da adduct.
2. `COMPOUND_STRUCTURE` gains `dimethyl_trisulfide` (class `trisulfide`, 2 carbons, not
   alpha,beta-unsaturated), which the module already names in `ADDUCT_POSITIVE_CLASSES` but has no
   structure record for.
3. The §2 header comment beside the k2 sec. B.3 35x aldehyde gap gains Sun 2025's **1118x**
   headspace-vs-fluorescence spread next to Bi 2022's 56x.
4. `data/keys/compounds.yml` gains `lenthionine` if the sulfur lane ever needs a cyclic polysulfide.
5. **No `REVERSIBLE_BINDING` row.** The sulfur gap stays open and is now *documented* rather than
   silently empty, with a measured magnitude bound of 1e-2 to 4e-2 L/g attached to it.

## 5. Flags

1. **Equation 1 is mislabelled and the label is dangerous.** "Flavor retention rate (%) =
   (H0i − H1i)/H0i × 100 %" with H0i = peak area **without** protein and H1i = **with** protein is
   the fraction of headspace signal that DISAPPEARED. That is a depletion/binding percentage. Bi
   2022 calls the identical quantity "binding percentage" and reserves "retention" for its PRV
   quantity. **A reader who takes Sun's "retention rate" as the fraction remaining in the headspace
   inverts every number in this paper.** The direction is unambiguous from the formula and from the
   paper's own reading of it ("PPI has adsorption effects"), so the values are safe; the name is not.
2. **There is no single binding temperature.** Equilibrium at **4 C for 12 h**, headspace sampled at
   **50 C for 40 min**. Those are different regimes: whatever bound in the cold has 40 minutes at
   50 C to re-partition during extraction, and any thermally activated covalent chemistry (which
   `matrix_sites.py` puts at Ea 15-23 kJ/mol) runs faster at 50 C than at 4 C. **The measured
   depletion is a hybrid of a cold equilibrium and a hot extraction, and the paper never separates
   them.**
3. **SPME depletion is not a partition measurement, and this is the single biggest obstacle to using
   the paper.** A CAR/DVB/PDMS fibre has finite, compound-dependent capacity; the three analytes
   compete with each other, with the **10 % v/v methanol** and with everything the isolate outgasses.
   Fibre uptake is proportional to gas concentration only in the non-depletive, non-saturated,
   non-competitive limit, which nobody demonstrates. **Every `K_g` in §3 item 3 rests on that
   unproven proportionality and every one of them is mine, not the paper's.** By contrast Bi 2022's
   shipped pea rows come from phase-ratio variation, an actual partition determination.
4. **There is no water leg printed as a number.** The control ("PBS solution in place of the PPI
   solution") is a matched leg and is exactly right in design, but its peak areas H0i live only in
   Fig. 1A's normalisation. The registry's ratio construction therefore has to be run backwards out
   of a percentage rather than forwards from two coefficients, which loses the error structure
   entirely: **no uncertainty on any reconstructed K_g can be stated.**
5. **10 % v/v methanol in every headspace and every zeta/CD/UV sample (mine).** Eight times Bi
   2022's ~1.25 %. Methanol competes for hydrophobic sites, changes the air/water partition of all
   three sulfur compounds, and at 10 % begins to perturb protein conformation in its own right. It
   cancels between Sun's sample and Sun's control. **It does not cancel against any other paper in
   the corpus**, so a Sun number and a Bi number are not on the same solvent basis. The CD blank is
   methanol, which shows the authors were aware of it there and nowhere else.
6. **KD = KS to the last decimal in all nine rows of Table 2, which cannot be right.** Equation 3 is
   `F0/F = (1 + KD[Q])(1 + KS[Q])`, a two-parameter fit whose whole purpose is to separate dynamic
   from static quenching. That the two parameters come out **exactly** equal — 319.71 and 319.71,
   744.70 and 744.70, 3558.45 and 3558.45, and so on — for three compounds at three temperatures is
   not a physical result. Either the polynomial was fitted with the constraint KD = KS (in which case
   it is a one-parameter fit and separates nothing, and `Ra`, a single correlation coefficient "of
   KD and KS", is consistent with that), or one column was duplicated in typesetting. **The paper's
   central claim in §3.4.2 — "static quenching as the primary mechanism" — is argued from the
   temperature dependence of KD and KS, and if KD ≡ KS by construction then that argument has no
   content.** Both columns should be treated as a single quenching constant KSV. **Ask the authors.**
7. **DMT's thermodynamic block does not close and its Ka column is not a van 't Hoff line.** Ka falls
   6.8x over the first six kelvin and 1.14x over the next six (mine); −RT ln(Ka) misses the printed
   ΔG by up to 1.47 kJ/mol where LEN misses by 0.13 (mine). ΔH = −131.74 kJ/mol is therefore set
   almost entirely by one point. **DMT's ΔH, ΔS and ΔG should not be transferred.** DMT is also the
   compound whose class (`trisulfide`) the repository already places in `ADDUCT_POSITIVE_CLASSES`,
   so an unmeasured irreversible channel is exactly what a collapsing "equilibrium constant" would
   look like.
8. **The Hill coefficient n is never printed, which leaves Ka uninterpretable.** Equation 4 fits
   `lg((F0−F)/F) = lg Ka + n lg[Q]`; Ka's units are L/mol only if n = 1, and its magnitude scales
   with whatever n was. LEN's 473 810 L/mol against DMD's 270 L/mol is a 1755x span that a modest
   difference in n would largely manufacture. **Request the n values.**
9. **Excitation at 260 nm is unusual and sits on the paper's own absorption peak.** The paper's UV
   section reports the PPI absorption maximum "around 260 nm" and reports that **it rises with SCFC
   concentration**. Exciting fluorescence at the wavelength whose absorbance is changing is exactly
   the inner-filter geometry that inflates apparent quenching. The authors did subtract SCFC-only
   blanks (which corrects the ligand's own absorbance and emission), but a **primary inner-filter
   correction for the changing sample absorbance at 260 nm is not described**. Standard practice
   (and Bi 2022's choice) is 290-295 nm.
10. **Nothing in this paper is a covalent measurement.** No mass spectrometry of the protein, no
    adduct search, no free-thiol assay before and after, no SDS-PAGE. The word "non-covalent" in the
    title is supported by three bond-breaker reagents read as significant/not-significant on a
    figure. **For DMT in particular, whose class the repository gates as adduct-positive, the paper
    supplies no control at all.** Compare the standing precedent: Amendment 6 ruling 4 sized 22-33 %
    of Meynier's t-2-hexenal "partition" as irreversible chemistry and quarantined the row.
11. **The bond-breaker run is at a different protein loading from the binding run.** 2 mg/mL against
    18 mg/mL — a **9x** difference (mine). The reversibility demonstration is therefore not made at
    the loading whose depletion is being explained.
12. **The protein is a black box.** Commercial isolate, "purity 90 %", supplier-stated. No thiol, no
    free amine, no SDS-PAGE, no solubility, no isoelectric point, no lot number.
    `protein_matrices.yml`'s `pea_isolate` site densities come from Gao 2020, Xiao 2024, Chen 2022,
    Shen 2022 and Chihi 2016 on *other* preparations, and Bi 2022's protein was made in-house.
    **Any pairing of a Sun constant with those site densities is a cross-preparation pairing.**
13. **The one odour threshold in the paper is a bare citation.** "0.008 ug/kg for dimethyl
    disulfide" (§1), attributed to J. Sun et al. 2021, with **no medium, no panel size, no
    orthonasal/retronasal distinction and no detection/recognition distinction**. It is not usable as
    a threshold record; chase the primary source.
14. **The E-nose is not sensory data.** §3.6 reports that SCFCs "amplify meat-like sulfur aroma and
    diminish bean-associated alcohol and aldehyde odors" on the basis of metal-oxide sensor response
    values. **There is no human panel anywhere in this paper.** W1W and W2W rising after adding
    sulfur compounds is a statement about sensor cross-sensitivity, not about perceived aroma, and
    the claim that beany notes were *reduced* rests on a W2S decrease that could equally be sensor
    competition. Nothing from §3.6 should enter any perceptual layer.
15. **The retention curve is non-monotone with a printed maximum at 0.5 mM and no explanation
    beyond site exposure.** Retention rises 0.1 -> 0.5 mM and falls above. The proposed mechanism —
    that rising flavour concentration unfolds the protein and *exposes* more sites — is the opposite
    of the usual saturation reading and is not tested. Bi 2022 saw the same shape and attributed the
    fall to exceeding the binding capacity. **The reconstructed K_g values in §3 item 3 are taken at
    the maximum, which is the most favourable point on the curve**; at 1.0 mM they would be lower and
    at 0.1 mM lower again. They are therefore an **upper** estimate within this paper's own range.
16. **The 0.5 mM values themselves are an inference, not a printed table.** §3.1 prints three ranges
    over 0.1-1.0 mM and states that all three peaked at 0.5 mM. The upper endpoint of each range is
    therefore the 0.5 mM value — a sound inference, but **the number is not printed against the
    concentration**, and the eight-point curve is Fig. 1A. Marked accordingly throughout.
17. **Lenthionine is named wrongly in the introduction.** The paper writes "1,2,3,5,6-pentathiolane";
    lenthionine is **1,2,3,5,6-pentathiepane** (a seven-membered ring, C2H4S5). The PubChem ID it
    gives (67,521) is lenthionine, so the compound is right and the name is a slip. Worth knowing
    before anyone keys it into `compounds.yml` from this paper's text.
18. **The GC column is non-polar (HP-5MS) where Bi 2022 used HP-WAX**, and full-scan 30-450 m/z where
    Bi used SIM. Peak-area precision on trace sulfur volatiles in full scan is worse than in SIM, and
    no limit of detection, limit of quantitation or replicate CV is reported anywhere.
19. **No error bars are printed on anything in Table 1 or Table 2.** The paper says triplicate and
    reports ANOVA letters on figures, but Table 1's H0 percentages and Table 2's nine Ka values carry
    **no standard deviation**. Bi 2022, by comparison, prints ±SE on every Klotz constant. **Every
    number in §4 above is a point with no stated dispersion.**
20. **What this paper does NOT contain**: any partition coefficient; any water-leg peak area as a
    number; any pH other than 7.2; any ionic-strength series beyond the single 0.6 M NaCl
    bond-breaker; any heat-treated protein; any covalent adduct measurement; any rate constant; any
    activation energy; any measured odour threshold; any human sensory panel; any thiol (so nothing
    directly on 2-methyl-3-furanthiol or 2-furfurylthiol); any Hill coefficient; any error bar on a
    tabulated constant; any protein molar mass.
21. **What to request from the authors**: (i) the eight-point retention table behind Fig. 1A with the
    control peak areas H0i, which would turn the reconstruction of §3 item 3 into a printed ratio;
    (ii) whether KD and KS were fitted independently, and if so why they are identical to two
    decimals in nine rows; (iii) the nine Hill coefficients n from Eq. 4; (iv) the bond-breaker peak
    areas of Fig. 4C as numbers, which would give the **size** of the reversible fraction for DMDS
    rather than only its sign; (v) whether any mass-spectrometric search for a DMDS or DMTS protein
    adduct was attempted; (vi) the SD on every cell of Tables 1 and 2; (vii) the Rq values and the
    zeta/particle-size values as numbers; (viii) whether the 12 h at 4 C equilibration was verified
    to reach equilibrium.
